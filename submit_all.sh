#!/usr/bin/env bash
# submit_all.sh
#
# Discovers all cam-ccd combos and orbits automatically, then submits
# the full pipeline for each, using reduced array sizes to stay within
# the 2000 job submit limit.
#
# Array sizes are set to N=6 (do_phot/copy_phot/cleanup_lc) and M=20
# (clean_phot), so that 16 cam-ccd x 8 orbits stays under 2000:
#   16 x (17 orbits x 2 x 6 + 6 + 20) = 1952 tasks
#
# Usage:
#   bash submit_all.sh <sector> <lcdir>
#
# Example:
#   bash submit_all.sh 100 lc
#   bash submit_all.sh 97 lc

set -eo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <sector> <lcdir>"
    exit 1
fi

sector=$1
lcdir=$2

sectoruse=$(printf "s%04d" "$sector")
dhome=$(pwd)
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
logdir="${LOG_DIR:-/lustre/research/mfausnau/logs}"

# ---- Sanity checks -----------------------------------------------------------
for var in DATA_DIR ISIS_DIR PIPELINE_DIR; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var is not set."
        exit 1
    fi
done
if [[ ! -d "$logdir" ]]; then
    echo "ERROR: Log directory does not exist: $logdir  (mkdir -p $logdir)"
    exit 1
fi

# ---- Discover cam-ccd combos from local phot.data files ----------------------
mapfile -t phot_files < <(ls "$dhome/sector${sector}/cam"*"_ccd"*/phot.data 2>/dev/null)

if [[ ${#phot_files[@]} -eq 0 ]]; then
    echo "ERROR: No sector${sector}/cam*_ccd*/phot.data files found in $dhome"
    exit 1
fi

echo "======================================================"
echo " full-sector pipeline submission"
echo " Sector  : $sectoruse"
echo " lcdir   : $lcdir"
echo " Found ${#phot_files[@]} cam-ccd combo(s)"
echo " Working directory: $dhome"
echo "======================================================"

# ---- Loop over each cam-ccd --------------------------------------------------
all_cleanup_lc_jobs=()
all_clean_phot_jobs=()
all_do_phot_job_ids=()
all_copy_phot_job_ids=()

for phot_file in "${phot_files[@]}"; do
    # Extract cam and ccd from path: sector100/cam2_ccd3/phot.data
    dir=$(dirname "$phot_file")
    camccd=$(basename "$dir")          # e.g. cam2_ccd3
    cam=${camccd#cam}                  # e.g. 2_ccd3
    cam=${cam%%_*}                     # e.g. 2
    ccd=${camccd#*_ccd}                # e.g. 3

    dtarget="${DATA_DIR}/${sectoruse}/cam${cam}-ccd${ccd}"

    if [[ ! -d "$dtarget" ]]; then
        echo "WARNING: Data directory not found, skipping: $dtarget"
        continue
    fi

    # Discover orbits that actually have slice directories
    mapfile -t orbits < <(
        for o_dir in "$dtarget"/o??; do
            [[ -d "$o_dir" ]] || continue
            if ls -d "$o_dir/slice"* &>/dev/null; then
                basename "$o_dir"
            fi
        done
    )

    if [[ ${#orbits[@]} -eq 0 ]]; then
        echo "WARNING: No orbits with slice dirs found for cam${cam}-ccd${ccd}, skipping."
        continue
    fi

    echo ""
    echo "------------------------------------------------------"
    echo " cam=$cam ccd=$ccd  orbits: ${orbits[*]}"
    echo "------------------------------------------------------"

    # -- Interactive cleanup (phot.data + aggregated lc dirs) --
    bash "$script_dir/cleanup_phot.sh" "$sector" "$cam" "$ccd" "$lcdir"

    # -- Submit cleanup_lc with reduced array size --
    cleanup_lc_job=$(sbatch \
        --parsable \
        --array=0-5%6 \
        --time=1-00:00:00 \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        "$script_dir/cleanup_lc.sbatch" \
        "$sector" "$cam" "$ccd" "$lcdir")
    echo "  cleanup_lc job ID: $cleanup_lc_job"
    all_cleanup_lc_jobs+=("$cleanup_lc_job")

    copy_job_ids=()

    for o in "${orbits[@]}"; do
        orbit_dir="$dtarget/$o"

        # Copy phot.data before any array task starts
        phot_dst="$orbit_dir/phot.data"
        echo "  Copying phot.data -> $phot_dst"
        cp "$phot_file" "$phot_dst"


        # Submit do_phot (reduced array, after cleanup_lc)
        do_job=$(sbatch \
            --parsable \
            --array=0-5%6 \
            --time=1-00:00:00 \
            --output="${logdir}/%x.o%A_%a" \
            --error="${logdir}/%x.e%A_%a" \
            --dependency=afterok:"$cleanup_lc_job" \
            "$script_dir/do_phot_em2.sbatch" \
            "$sector" "$cam" "$ccd" "$o" "$lcdir")
        echo "    do_phot  [$o] job ID: $do_job"

        # Submit copy_phot (reduced array, after do_phot)
        copy_job=$(sbatch \
            --parsable \
            --array=0-5%6 \
            --time=1-00:00:00 \
            --output="${logdir}/%x.o%A_%a" \
            --error="${logdir}/%x.e%A_%a" \
            --dependency=afterany:"$do_job" \
            "$script_dir/copy_phot_em2.sbatch" \
            "$sector" "$cam" "$ccd" "$o" "$lcdir")
        echo "    copy_phot[$o] job ID: $copy_job"

        copy_job_ids+=("$copy_job")
        all_do_phot_job_ids+=("$do_job")
        all_copy_phot_job_ids+=("$copy_job")
    done

    # -- Submit clean_phot after all copy_phot for this cam-ccd --
    dep_str="afterany:$(IFS=:; echo "${copy_job_ids[*]}")"
    clean_job=$(sbatch \
        --parsable \
        --array=1-20%20 \
        --time=1-00:00:00 \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        --dependency="$dep_str" \
        "$script_dir/clean_phot_em2.sbatch" \
        "$sector" "$cam" "$ccd" "$lcdir")
    echo "  clean_phot job ID: $clean_job"
    all_clean_phot_jobs+=("$clean_job")

done

# ---- Submit single aggregate_errors job after all clean_phot jobs -----------
echo ""
echo "Submitting aggregate_errors (after all clean_phot jobs)..."
all_clean_dep="afterany:$(IFS=:; echo "${all_clean_phot_jobs[*]}")"
timestamp=$(date +%Y%m%d_%H%M%S)
agg_job=$(sbatch \
    --parsable \
    --time=1-00:00:00 \
    --output="${logdir}/%x.o%A" \
    --error="${logdir}/%x.e%A" \
    --dependency="$all_clean_dep" \
    "$script_dir/aggregate_errors.sbatch" \
    "$sector" "all" "all" "$lcdir" "${all_do_phot_job_ids[@]}" "--" "${all_copy_phot_job_ids[@]}" "--" "${all_clean_phot_jobs[@]}")
echo "  aggregate_errors job ID: $agg_job"

# ---- Submit purge_lc after aggregate_errors ----------------------------------
echo ""
echo "Submitting purge_lc (after aggregate_errors)..."
purge_job=$(sbatch \
    --parsable \
    --time=1-00:00:00 \
    --output="${logdir}/%x.o%A_%a" \
    --error="${logdir}/%x.e%A_%a" \
    --dependency=afterany:"$agg_job" \
    "$script_dir/purge_lc.sbatch" \
    "$sector" "all" "all" "$lcdir")
echo "  purge_lc job ID: $purge_job"

# ---- Summary -----------------------------------------------------------------
echo ""
echo "======================================================"
echo " All cam-ccd combos submitted."
echo " cleanup_lc jobs : ${all_cleanup_lc_jobs[*]}"
echo " clean_phot jobs : ${all_clean_phot_jobs[*]}"
echo " purge_lc job    : $purge_job"
echo ""
echo " Monitor with:"
echo "   squeue -u \$USER"
echo "   ls $logdir/{do_phot,copy_phot,clean_phot,cleanup_lc}.o*"
echo "   cat errorlog-$(printf 's%04d' $sector)-camall-ccdall-*.txt  (after pipeline completes)"
echo "======================================================"
