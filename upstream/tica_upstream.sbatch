#!/usr/bin/env bash
#SBATCH --job-name=tica_upstream
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=2-00:00:00
#SBATCH --array=0-15%4

# HPCC upstream processing (download -> subtract -> slice tarballs) for one
# sector, one CCD per array task: task i -> cam i/4+1, ccd i%4+1, so
#   0-3 cam1 ccd1-4, 4-7 cam2, 8-11 cam3, 12-15 cam4   (cam3 ccd3 = task 10)
# Each task runs upstream/tica_upstream.sh with one worker per CPU.
#
# Usage (after: source setup_vars_hpcc.sh):
#   sbatch $TESS_CODE_DIR/upstream/tica_upstream.sbatch <sector> [orbit|all]
#   sbatch --array=10 $TESS_CODE_DIR/upstream/tica_upstream.sbatch 102 o1a   # one CCD
#
# %4 keeps at most 4 CCDs (64 concurrent MAST downloads) running at once; raise
# it with care. If compute nodes have no internet, prefetch on a login node
# first: python $TESS_CODE_DIR/upstream/tica_prefetch.py ... (see that script).

set -eo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch [--array=...] $0 <sector> [orbit|all]"
    exit 1
fi
sector=$1
orbit=${2:-all}
task=${SLURM_ARRAY_TASK_ID:-0}
cam=$(( task / 4 + 1 ))
ccd=$(( task % 4 + 1 ))

if [[ -z "${TESS_CODE_DIR:-}" ]]; then
    echo "ERROR: TESS_CODE_DIR is not set; run 'source setup_vars_hpcc.sh' before sbatch"
    exit 1
fi
source "${SETUP_VARS:-$TESS_CODE_DIR/setup_vars_hpcc.sh}"
export NPROC=${SLURM_CPUS_PER_TASK:-$NPROC}

echo "tica_upstream: sector $sector cam $cam ccd $ccd orbit $orbit on $(hostname), $NPROC workers"
bash "$TESS_CODE_DIR/upstream/tica_upstream.sh" "$sector" "$cam" "$ccd" "$orbit"
