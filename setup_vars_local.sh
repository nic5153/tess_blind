#!/usr/bin/env bash
# Local (WSL) replacement for ~/setup_vars.sh on HPCC. Source it before running
# anything:   source setup_vars_local.sh
#
#   DATA_DIR      upstream products (the HPCC ISIS data directory's role)
#   WORK_DIR      role of /home/nimcclur/work/TESS/photometry (sector*/ outputs);
#                 run orbit_pipeline.sbatch from here
#   TESS_CODE_DIR this repository (scripts are called from here, not copied)
#
# It also puts local_slurm/ (sbatch/squeue stand-ins) first on PATH and selects
# the Python photometry (PHOT_BACKEND=python) in place of ISIS phot2.

TESS_CODE_DIR=${TESS_CODE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}
export TESS_CODE_DIR
export SETUP_VARS=${SETUP_VARS:-$TESS_CODE_DIR/setup_vars_local.sh}
export TESS_LOCAL_ROOT=${TESS_LOCAL_ROOT:-/mnt/d/tess_blind_data}

# Separate trees per compression setting, so lossless and lossy runs can be
# compared side by side (tica/ + photometry/ vs tica_q4/ + photometry_q4/).
suffix=${COMPRESS_Q:+_q${COMPRESS_Q}}
export DATA_DIR=${DATA_DIR:-$TESS_LOCAL_ROOT/tica$suffix}
export WORK_DIR=${WORK_DIR:-$TESS_LOCAL_ROOT/photometry$suffix}
unset suffix
export LOG_DIR=${LOG_DIR:-$WORK_DIR/logs}
export RAW_CACHE=${RAW_CACHE:-$TESS_LOCAL_ROOT/raw_cache}
export TICA_URL=${TICA_URL:-https://archive.stsci.edu/hlsps/tica}

# Downstream: Python photometry instead of ISIS, local SLURM stand-ins, and an
# optional tesstools/ folder (analysis.sbatch skips matchbkgtolc.py without it).
export PHOT_BACKEND=${PHOT_BACKEND:-python}
export LOCAL_SLURM=1
# Delete each orbit's difference-image tarballs once orbit_pipeline has
# finished it. Set PURGE_IMAGES=0 to keep them (e.g. to rerun photometry or
# compare compression settings on the images themselves).
export PURGE_IMAGES=${PURGE_IMAGES:-1}
export TESSTOOLS_DIR=${TESSTOOLS_DIR:-$TESS_CODE_DIR/tesstools}
case ":$PATH:" in
    *":$TESS_CODE_DIR/local_slurm:"*) ;;
    *) export PATH="$TESS_CODE_DIR/local_slurm:$PATH" ;;
esac

# Parallel workers (upstream slices and local array tasks). Each needs ~0.5 GB;
# WSL currently has ~7 GB.
export NPROC=${NPROC:-4}
# copy_phot tasks all append (cat >>) to the same aggregate lc files; appends on
# the WSL D: mount are not guaranteed atomic, so run that step as one task.
export COPY_PHOT_ARRAY=${COPY_PHOT_ARRAY:-0-0}
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1

# ---- upstream processing knobs (defaults shown) -------------------------------
# SLICE_SIZE=auto        frames per slice; auto = round(N_halforbit/450) like HPCC
# HALF_SPLIT=midpoint    split orbits without a mid-orbit gap at the midpoint
# REF_NFRAMES=20         lowest-background frames median-combined into ref.fits
# REF_NCAND=60           candidates screened for the reference
# SUB_MODE=reference     reference | adjacent (frame minus its neighbor in time)
# KERNEL_METHOD=lstsq    lstsq (real-space, ISIS-like) | photutils (Fourier matching)
# KERNEL_HALF=2          lstsq kernel half-width: 2 -> 5x5 kernel
# NTILE=4                NTILE x NTILE kernel/PSF tiles across the CCD
# KSIZE=13               PSF / kernel stamp size [pix]
# DIFF_BKG_BOX=16        Background2D box for sky removal on the difference image
# COLUMN_CLEAN=1         remove column offsets, zero noisy (bleed/smear) columns
# COLMASK_NSIG=5         a column is zeroed if its noise > this x the typical column
# COLMASK_GROW=1         also zero this many neighboring columns on each side
# SKY_BOX=64             block size of the frame sky map stored in bkg_<name>
# INTERP_MODE=header     header | full: interp_ carries only the WCS header, or the
#                        full registered frame (doubles disk use)
# COMPRESS_Q=            fpack RICE quantize level for lossy test (empty = off)
# RMS_QUAL_MASK=-1       QUAL_BIT bits that exclude a frame from rms_list_tica
# KEEP_RAW=0             1 keeps downloaded FFIs (compressed .fz when COMPRESS_Q set)
# EPSF_OVERSAMP=2        ePSF oversampling for the PSF photometry
# ---- photometry knobs (phot_slice.py) -----------------------------------------
# PHOT_FIT_SHAPE=5       PSF fit box [pix]
# PHOT_XY_OFFSET=0.5     added to phot.data positions (do_photutils_extract.py
#                        stores centroid - 0.5)

if [[ -f "$HOME/miniconda3/etc/profile.d/conda.sh" && "${CONDA_DEFAULT_ENV:-}" != "tess" ]]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
    conda activate tess
fi

mkdir -p "$DATA_DIR" "$WORK_DIR" "$LOG_DIR" "$RAW_CACHE"
