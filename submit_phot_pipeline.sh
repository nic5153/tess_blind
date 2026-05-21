#!/usr/bin/env bash
# submit_phot_pipeline.sh
#
# Submits the full photometry pipeline for one cam-ccd combination,
# across one or more orbit directories.  Handles dependency chaining:
#
#   cleanup_phot (interactive: phot.data + aggregated lc dirs)
#   cleanup_lc   (sbatch array: per-slice lc dirs)
#      |
#      +-- for each orbit: do_phot -> copy_phot
#                                          |
#                                     clean_phot (after ALL copy_phot finish)
#
# Usage:
#   bash submit_phot_pipeline.sh <sector> <cam> <ccd> <lcdir> <orbit1> [orbit2 ...]
#
# Examples:
#   # Typical sector, two orbits:
#   bash submit_phot_pipeline.sh 100 1 1 lc o1a o1b
#
#   # Special sector, all 8 orbits:
#   bash submit_phot_pipeline.sh 97 2 3 lc o1a o1b o2a o2b o3a o3b o4a o4b
#
# Notes:
#   - Must be run from the phot_scripts directory (where sector*/cam*_ccd* live).
#   - DATA_DIR, ISIS_DIR, and PIPELINE_DIR must be set in your environment.
#   - Logs go to /lustre/research/mfausnau/logs/ (must exist).

set -eo pipefail

# ---- Argument parsing --------------------------------------------------------
if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <sector> <cam> <ccd> <lcdir> <orbit1> [orbit2 ...]"
    echo ""
    echo "  sector  : integer sector number (e.g. 100)"
    echo "  cam     : camera number (1-4)"
    echo "  ccd     : ccd number (1-4)"
    echo "  lcdir   : light curve subdirectory name (e.g. lc)"
    echo "  orbitN  : one or more orbit labels (e.g. o1a o1b o2a o2b)"
    exit 1
fi

sector=$1
cam=$2
ccd=$3
lcdir=$4
shift 4
orbits=("$@")

sectoruse=$(printf "s%04d" "$sector")
dhome=$(pwd)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir="${LOG_DIR:-/lustre/research/mfausnau/logs}"

# ---- Sanity checks -----------------------------------------------------------
if [[ -z "${DATA_DIR:-}" ]]; then
    echo "ERROR: DATA_DIR is not set."
    exit 1
fi
if [[ -z "${ISIS_DIR:-}" ]]; then
    echo "ERROR: ISIS_DIR is not set."
    exit 1
fi
if [[ -z "${PIPELINE_DIR:-}" ]]; then
    echo "ERROR: PIPELINE_DIR is not set."
    exit 1
fi
if [[ ! -d "$logdir" ]]; then
    echo "ERROR: Log directory does not exist: $logdir"
    echo "       Create it with: mkdir -p $logdir"
    exit 1
fi

dtarget="${DATA_DIR}/${sectoruse}/cam${cam}-ccd${ccd}"
if [[ ! -d "$dtarget" ]]; then
    echo "ERROR: Data directory not found: $dtarget"
    exit 1
fi

echo "======================================================"
echo " Photometry pipeline submission"
echo " Sector : $sectoruse   Cam: $cam   CCD: $ccd"
echo " lcdir  : $lcdir"
echo " Orbits : ${orbits[*]}"
echo " Working directory: $dhome"
echo "======================================================"

# ---- Step 1: Interactive pre-flight cleanup (phot.data + aggregated lc dirs) ---
echo ""
echo "[1/5] Running interactive pre-flight cleanup (phot.data + aggregated lc dirs)..."
bash "$script_dir/cleanup_phot.sh" "$sector" "$cam" "$ccd" "$lcdir"

# ---- Step 2: Submit cleanup_lc (per-slice lc dir removal) -------------------
echo ""
echo "[2/5] Submitting cleanup_lc array job..."
cleanup_lc_job=$(sbatch \
    --parsable \
    --output="${logdir}/%x.o%A_%a" \
    --error="${logdir}/%x.e%A_%a" \
    "$script_dir/cleanup_lc.sbatch" \
    "$sector" "$cam" "$ccd" "$lcdir")
echo "  cleanup_lc job ID: $cleanup_lc_job"

# ---- Step 3: Submit do_phot + copy_phot for each orbit ----------------------
echo ""
echo "[3/5] Submitting do_phot -> copy_phot chains per orbit (after cleanup_lc)..."

copy_job_ids=()    # collect all copy_phot job IDs for the final dependency
do_phot_job_ids=()  # collect all do_phot job IDs for error aggregation

for o in "${orbits[@]}"; do

    # Verify the orbit directory exists
    orbit_dir="$dtarget/$o"
    if [[ ! -d "$orbit_dir" ]]; then
        echo "  WARNING: Orbit directory not found, skipping: $orbit_dir"
        continue
    fi

    n_slices=$(ls -d "$orbit_dir/slice"* 2>/dev/null | wc -l)
    echo ""
    echo "  Orbit: $o  ($n_slices slices)"

    # Copy phot.data into place before any array task starts
    phot_src="$dhome/sector${sector}/cam${cam}_ccd${ccd}/phot.data"
    phot_dst="$orbit_dir/phot.data"
    if [[ ! -e "$phot_src" ]]; then
        echo "  ERROR: phot.data source not found: $phot_src"
        exit 1
    fi
    echo "  Copying phot.data -> $phot_dst"
    cp "$phot_src" "$phot_dst"


    # Submit do_phot, dependent on cleanup_lc completing successfully
    do_job=$(sbatch \
        --parsable \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        --dependency=afterok:"$cleanup_lc_job" \
        "$script_dir/do_phot_em2.sbatch" \
        "$sector" "$cam" "$ccd" "$o" "$lcdir")
    echo "    do_phot  job ID: $do_job  (after cleanup_lc $cleanup_lc_job)"

    # Submit copy_phot, dependent on do_phot completing successfully
    copy_job=$(sbatch \
        --parsable \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        --dependency=afterany:"$do_job" \
        "$script_dir/copy_phot_em2.sbatch" \
        "$sector" "$cam" "$ccd" "$o" "$lcdir")
    echo "    copy_phot job ID: $copy_job  (after do_phot $do_job)"

    copy_job_ids+=("$copy_job")
    do_phot_job_ids+=("$do_job")
done

if [[ ${#copy_job_ids[@]} -eq 0 ]]; then
    echo "ERROR: No orbits were submitted (no valid orbit directories found)."
    exit 1
fi

# ---- Step 3: Submit clean_phot after ALL copy_phot jobs finish ---------------
echo ""
echo "[4/5] Submitting clean_phot (after all copy_phot jobs)..."

# Build colon-separated dependency string: afterany:id1:id2:id3...
dep_str="afterany:$(IFS=:; echo "${copy_job_ids[*]}")"

clean_job=$(sbatch \
    --parsable \
    --output="${logdir}/%x.o%A_%a" \
    --error="${logdir}/%x.e%A_%a" \
    --dependency="$dep_str" \
    "$script_dir/clean_phot_em2.sbatch" \
    "$sector" "$cam" "$ccd" "$lcdir")
echo "    clean_phot job ID: $clean_job  (after ${copy_job_ids[*]})"

# ---- Step 5: Submit aggregate_errors after clean_phot ------------------------
echo ""
echo "[5/6] Submitting aggregate_errors (after clean_phot)..."
agg_job=$(sbatch \
    --parsable \
    --output="${logdir}/%x.o%A" \
    --error="${logdir}/%x.e%A" \
    --dependency=afterany:"$clean_job" \
    "$script_dir/aggregate_errors.sbatch" \
    "$sector" "$cam" "$ccd" "$lcdir" "${do_phot_job_ids[@]}" "--" "${copy_job_ids[@]}" "--" "$clean_job")
echo "    aggregate_errors job ID: $agg_job"

# ---- Step 6: Submit purge_lc after aggregate_errors --------------------------
echo ""
echo "[6/6] Submitting purge_lc (after aggregate_errors)..."
purge_job=$(sbatch \
    --parsable \
    --output="${logdir}/%x.o%A_%a" \
    --error="${logdir}/%x.e%A_%a" \
    --dependency=afterany:"$agg_job" \
    "$script_dir/purge_lc.sbatch" \
    "$sector" "$cam" "$ccd" "$lcdir")
echo "    purge_lc job ID: $purge_job"

# ---- Summary -----------------------------------------------------------------
echo ""
echo "======================================================"
echo " Pipeline submitted successfully."
echo ""
echo " Job chain:"
echo "   cleanup_lc[$cleanup_lc_job]"
for i in "${!orbits[@]}"; do
    echo "     -> ${orbits[$i]}: do_phot -> copy_phot[${copy_job_ids[$i]}]"
done
echo "   All copy_phot -> clean_phot[$clean_job] -> aggregate_errors[$agg_job] -> purge_lc[$purge_job]"
echo ""
echo " Monitor with:"
echo "   squeue -u \$USER"
echo "   # Follow a specific array task log (task 1 of first do_phot job):"
echo "   tail -f $logdir/do_phot.o${do_job}_1"
echo "   # List all logs for this pipeline run:"
echo "   ls $logdir/{do_phot,copy_phot,clean_phot,cleanup_lc}.o*"
echo "======================================================"
