#!/usr/bin/env bash
# HPCC (TTU nocona) environment for the self-processed TICA pipeline: you
# download and subtract the FFIs yourself (upstream/) instead of using someone
# else's ISIS products. Source it on a login node
# before submitting anything; sbatch passes the environment on to the jobs.
#
#   source setup_vars_hpcc.sh
#
# Everything lives in one folder you own, laid out like the local version
# (D:\tess_blind_data on the PC):
#
#   TESS_ROOT  = /lustre/work/$USER/tess_blind_data
#   DATA_DIR   = $TESS_ROOT/tica        upstream products (ref.fits, slice
#                                       tarballs, rms_list_tica)
#   RAW_CACHE  = $TESS_ROOT/raw_cache   downloaded FFIs, deleted as each slice
#                                       is processed
#   WORK_DIR   = $TESS_ROOT/photometry  where orbit_pipeline runs (sector*/)
#   LOG_DIR    = $WORK_DIR/logs         SLURM logs
#
# /lustre/work is not purged but has a quota, so each half-orbit's difference
# images are deleted once orbit_pipeline has analysed it (PURGE_IMAGES=1), and
# the upstream disk guard checks the quota. See your quota with:
#   lfs quota -h -u $USER /lustre/work
# The code itself stays wherever you keep it (TESS_CODE_DIR, this folder).

export TESS_CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SETUP_VARS="$TESS_CODE_DIR/setup_vars_hpcc.sh"
export TESS_ROOT=${TESS_ROOT:-/lustre/work/$USER/tess_blind_data}

# All paths are derived from TESS_ROOT and always overwritten here: a shell that
# already loaded the old ~/setup_vars.sh has DATA_DIR pointing at the shared
# ISIS data, which these scripts must never write into. To relocate everything,
# set TESS_ROOT before sourcing this file.
# Separate trees per compression setting (tica/ + photometry/ vs tica_q4/ + photometry_q4/).
suffix=${COMPRESS_Q:+_q${COMPRESS_Q}}
export DATA_DIR="$TESS_ROOT/tica$suffix"
export WORK_DIR="$TESS_ROOT/photometry$suffix"
unset suffix
export RAW_CACHE="$TESS_ROOT/raw_cache"
export LOG_DIR="$WORK_DIR/logs"
export TICA_URL=${TICA_URL:-https://archive.stsci.edu/hlsps/tica}

# Python photometry instead of ISIS phot2; real SLURM (no local stand-ins).
export PHOT_BACKEND=python
export LOCAL_SLURM=0
# Delete each orbit's difference images once it is analysed (work quota).
# Set PURGE_IMAGES=0 to keep them.
export PURGE_IMAGES=${PURGE_IMAGES:-1}

# Workers per job: the job's CPU count inside SLURM, else a modest default.
export NPROC=${SLURM_CPUS_PER_TASK:-${NPROC:-8}}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# Python environment: your existing conda env "nmcclur" (needs photutils >= 2.0,
# astropy, scipy, astroquery, tess-point). Set TESS_CONDA_ENV to use another one.
export TESS_CONDA_ENV=${TESS_CONDA_ENV:-nmcclur}
if [[ "${CONDA_DEFAULT_ENV:-}" != "$TESS_CONDA_ENV" ]] && command -v conda > /dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate "$TESS_CONDA_ENV" || echo "WARNING: conda env $TESS_CONDA_ENV not found (conda info --envs)" >&2
fi

mkdir -p "$DATA_DIR" "$RAW_CACHE" "$WORK_DIR" "$LOG_DIR"
