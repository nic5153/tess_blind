#!/usr/bin/env bash
# HPCC (TTU nocona) environment for the self-processed TICA pipeline: you
# download and subtract the FFIs yourself (upstream/) instead of using someone
# else's ISIS products. Source it on a login node
# before submitting anything; sbatch passes the environment on to the jobs.
#
#   source setup_vars_hpcc.sh
#
#   DATA_DIR   your upstream products (ref.fits, slice tarballs, rms_list_tica)
#   RAW_CACHE  downloaded FFIs, deleted as each slice is processed
#   WORK_DIR   where orbit_pipeline runs (sector*/ outputs), as before
#   LOG_DIR    SLURM logs, as before
#
# DATA_DIR and RAW_CACHE default to /lustre/scratch (large; files left unused
# for a long time are purged by HPCC policy). Check quotas with:
#   lfs quota -h -u $USER /lustre/scratch ; lfs quota -h -u $USER /lustre/work
# To use the old ISIS products instead, source ~/setup_vars.sh as before.

TESS_CODE_DIR=${TESS_CODE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
export TESS_CODE_DIR
export SETUP_VARS=${SETUP_VARS:-$TESS_CODE_DIR/setup_vars_hpcc.sh}

# Separate trees per compression setting (tica/ vs tica_q4/, ...).
suffix=${COMPRESS_Q:+_q${COMPRESS_Q}}
export DATA_DIR=${DATA_DIR:-/lustre/scratch/$USER/TESS/tica$suffix}
export RAW_CACHE=${RAW_CACHE:-/lustre/scratch/$USER/TESS/raw_cache}
export WORK_DIR=${WORK_DIR:-/home/$USER/work/TESS/photometry$suffix}
unset suffix
export LOG_DIR=${LOG_DIR:-/lustre/work/$USER/TESS/photometry/logs}
export TICA_URL=${TICA_URL:-https://archive.stsci.edu/hlsps/tica}

# Python photometry instead of ISIS phot2; real SLURM (no local stand-ins).
export PHOT_BACKEND=${PHOT_BACKEND:-python}
export LOCAL_SLURM=0
# Keep difference images after an orbit is analysed unless asked (HPCC has room).
export PURGE_IMAGES=${PURGE_IMAGES:-0}

# Workers per job: the job's CPU count inside SLURM, else a modest default.
export NPROC=${SLURM_CPUS_PER_TASK:-${NPROC:-8}}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Python environment (see flow_readme.txt for creating it once).
if [[ "${CONDA_DEFAULT_ENV:-}" != "${TESS_CONDA_ENV:-tess}" ]] && command -v conda > /dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "${TESS_CONDA_ENV:-tess}" || echo "WARNING: conda env ${TESS_CONDA_ENV:-tess} not found" >&2
fi

mkdir -p "$DATA_DIR" "$RAW_CACHE" "$WORK_DIR" "$LOG_DIR"
