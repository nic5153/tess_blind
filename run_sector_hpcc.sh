#!/usr/bin/env bash
# Submit a sector on HPCC as a chain of SLURM jobs, one half-orbit at a time:
#   upstream (one job per CCD)  ->  orbit_pipeline  ->  next half-orbit ...
# Run on a login node after:  source setup_vars_hpcc.sh
#
# Usage:
#   bash run_sector_hpcc.sh --sector N [--cam C] [--ccd D] [--orbit O] [--parallel P] [--dry-run]
#
#   --sector N    TESS sector (required)
#   --cam C       camera 1-4, or a comma list like 1,3   (default: all cameras)
#   --ccd D       CCD 1-4, or a comma list like 2,3      (default: all CCDs)
#   --orbit O     half-orbit(s) o1a,o1b,o2a,o2b          (default: all four)
#   --parallel P  CCDs processed at the same time         (default: 4)
#   --dry-run     print the sbatch commands without submitting
#
# Examples:
#   bash run_sector_hpcc.sh --sector 102                    # whole sector
#   bash run_sector_hpcc.sh --sector 102 --cam 3 --ccd 3    # one CCD
#   bash run_sector_hpcc.sh --sector 102 --cam 3            # all CCDs of camera 3
#
# Each orbit_pipeline waits for its upstream jobs (afterany: CCDs whose upstream
# failed are skipped or left unpurged), and each half-orbit's upstream waits for
# the previous orbit_pipeline, so with PURGE_IMAGES=1 only one half-orbit of
# difference images exists at a time. Rerunning resumes: finished slices are
# skipped. Each CCD needs ~50 GB per half-orbit at 200-s cadence.

set -eo pipefail

usage() {
    sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
}

sector="" cams="1,2,3,4" ccds="1,2,3,4" orbits="o1a,o1b,o2a,o2b" parallel=4 dry_run=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sector)     sector=$2; shift 2 ;;
        --sector=*)   sector=${1#*=}; shift ;;
        --cam)        cams=$2; shift 2 ;;
        --cam=*)      cams=${1#*=}; shift ;;
        --ccd)        ccds=$2; shift 2 ;;
        --ccd=*)      ccds=${1#*=}; shift ;;
        --orbit)      orbits=$2; shift 2 ;;
        --orbit=*)    orbits=${1#*=}; shift ;;
        --parallel)   parallel=$2; shift 2 ;;
        --parallel=*) parallel=${1#*=}; shift ;;
        --dry-run)    dry_run=1; shift ;;
        -h|--help)    usage; exit 0 ;;
        *) echo "ERROR: unknown option '$1'"; echo; usage; exit 1 ;;
    esac
done

if [[ ! "$sector" =~ ^[0-9]+$ ]]; then
    echo "ERROR: --sector is required (a number)"; echo; usage; exit 1
fi
if [[ ! "$parallel" =~ ^[0-9]+$ || "$parallel" -lt 1 ]]; then
    echo "ERROR: --parallel must be a positive number"; exit 1
fi

# cam/ccd pairs -> upstream array tasks (task = 4*(cam-1) + ccd-1, decoded by
# upstream/tica_upstream.sbatch).
tasks=() labels=()
for cam in ${cams//,/ }; do
    for ccd in ${ccds//,/ }; do
        if [[ ! "$cam" =~ ^[1-4]$ || ! "$ccd" =~ ^[1-4]$ ]]; then
            echo "ERROR: --cam and --ccd take values 1-4 (got cam=$cam ccd=$ccd)"; exit 1
        fi
        tasks+=($(( 4 * (cam - 1) + ccd - 1 )))
        labels+=("cam${cam}-ccd${ccd}")
    done
done
array="$(IFS=,; echo "${tasks[*]}")%${parallel}"

for orbit in ${orbits//,/ }; do
    if [[ ! "$orbit" =~ ^o[0-9][ab]$ ]]; then
        echo "ERROR: --orbit values look like o1a, o1b, o2a, o2b (got $orbit)"; exit 1
    fi
done

if [[ -z "${TESS_CODE_DIR:-}" || -z "${TESS_ROOT:-}" ]]; then
    echo "ERROR: run 'source setup_vars_hpcc.sh' first"
    exit 1
fi
# Refuse to submit if the data or results paths point outside your own
# TESS_ROOT (e.g. a shell where the old ~/setup_vars.sh was sourced afterwards).
for path in "$DATA_DIR" "$WORK_DIR" "$RAW_CACHE"; do
    if [[ "$path" != "$TESS_ROOT"/* ]]; then
        echo "ERROR: $path is outside TESS_ROOT=$TESS_ROOT."
        echo "       Run 'source setup_vars_hpcc.sh' again (after any other setup file) and retry."
        exit 1
    fi
done

submit() {  # print in dry-run mode, otherwise run sbatch and return the job id
    if [[ $dry_run -eq 1 ]]; then
        echo "sbatch $*" >&2
        echo "DRYRUN$RANDOM"
    else
        sbatch "$@"
    fi
}

echo "Sector $sector: ${#tasks[@]} CCD(s) [${labels[*]}], half-orbits ${orbits//,/ }, $parallel at a time"
echo "  DATA_DIR=$DATA_DIR"
echo "  WORK_DIR=$WORK_DIR  PURGE_IMAGES=$PURGE_IMAGES"
[[ $dry_run -eq 1 ]] || mkdir -p "$LOG_DIR" "$WORK_DIR"

prev=""
for orbit in ${orbits//,/ }; do
    up=$(submit --parsable --array="$array" ${prev:+--dependency=afterany:$prev} \
        --output="$LOG_DIR/%x.o%A_%a" --error="$LOG_DIR/%x.e%A_%a" \
        "$TESS_CODE_DIR/upstream/tica_upstream.sbatch" "$sector" "$orbit")
    pipe=$(submit --parsable --dependency=afterany:"$up" --chdir="$WORK_DIR" \
        --output="$LOG_DIR/%x.o%j" --error="$LOG_DIR/%x.e%j" \
        "$TESS_CODE_DIR/orbit_pipeline.sbatch" "$sector" "$orbit")
    echo "  $orbit: upstream job $up -> orbit_pipeline job $pipe"
    prev=$pipe
done
[[ $dry_run -eq 1 ]] && echo "(dry run: nothing was submitted)"
echo "Monitor with: squeue -u $USER   (logs in $LOG_DIR)"
