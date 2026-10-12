#!/usr/bin/env bash
# Submit a whole sector on HPCC as a chain of SLURM jobs, one half-orbit at a
# time:  upstream (array of CCDs)  ->  orbit_pipeline  ->  next half-orbit ...
# Run on a login node after:  source setup_vars_hpcc.sh
#
# Usage:
#   bash run_sector_hpcc.sh <sector> [ccd-array]
# Examples:
#   bash run_sector_hpcc.sh 102            # all 16 CCDs (array 0-15%4)
#   bash run_sector_hpcc.sh 102 10         # cam3 ccd3 only (task = 4*(cam-1) + ccd-1)
#   bash run_sector_hpcc.sh 102 8-11%2     # camera 3
#
# Each orbit_pipeline waits for its upstream jobs (afterany: it skips CCDs whose
# upstream failed), and each half-orbit's upstream waits for the previous
# orbit_pipeline, so with PURGE_IMAGES=1 only one half-orbit of difference
# images exists at a time. Rerunning resumes: finished slices are skipped.

set -eo pipefail

if [[ $# -lt 1 || $# -gt 2 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: bash $0 <sector> [ccd-array, default 0-15%4]"
    exit 1
fi
sector=$1
array=${2:-0-15%4}

if [[ -z "${TESS_CODE_DIR:-}" ]]; then
    echo "ERROR: run 'source setup_vars_hpcc.sh' first"
    exit 1
fi
mkdir -p "$LOG_DIR" "$WORK_DIR"

echo "Submitting sector $sector, CCD array $array"
echo "  DATA_DIR=$DATA_DIR  WORK_DIR=$WORK_DIR  PURGE_IMAGES=$PURGE_IMAGES"
prev=""
for orbit in ${ORBITS:-o1a o1b o2a o2b}; do
    up=$(sbatch --parsable --array="$array" ${prev:+--dependency=afterany:$prev} \
        --output="$LOG_DIR/%x.o%A_%a" --error="$LOG_DIR/%x.e%A_%a" \
        "$TESS_CODE_DIR/upstream/tica_upstream.sbatch" "$sector" "$orbit")
    pipe=$(sbatch --parsable --dependency=afterany:"$up" --chdir="$WORK_DIR" \
        --output="$LOG_DIR/%x.o%j" --error="$LOG_DIR/%x.e%j" \
        "$TESS_CODE_DIR/orbit_pipeline.sbatch" "$sector" "$orbit")
    echo "  $orbit: upstream job $up -> orbit_pipeline job $pipe"
    prev=$pipe
done
echo "Monitor with: squeue -u $USER   (logs in $LOG_DIR)"
