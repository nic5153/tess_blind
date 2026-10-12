FLOW FOR RMS SOURCE AND PHOTOMETRY

HPCC CODE UPDATE SAFETY:

The HPCC photometry directory also contains run products, logs, plot tarballs,
and recovered review files. Avoid switching branches there unless you know the
working tree is clean and disposable.

Safer update path for HPCC:

bash sync_code_from_codex.sh

This fetches GitHub and restores only code/script files from
origin/updated-hpcc (pass another ref to override), leaving sector outputs and
plot exports in place.

HPCC WITH YOUR OWN TICA PROCESSING (branch updated-hpcc):

Download and subtract the FFIs yourself (same code as the local version) under
SLURM, instead of relying on someone else's ISIS products.

One-time Python environment on HPCC (in a login shell):
conda create -n tess -c conda-forge python=3.12 astropy photutils astroquery scipy numpy pandas matplotlib tqdm
conda activate tess && pip install tess-point

Check whether compute nodes can reach MAST (prints an HTTP status line if so):
srun -p nocona -t 00:05:00 curl -sI https://archive.stsci.edu/hlsps/tica/ | head -1
If not, prefetch each CCD half-orbit on a login node before its upstream job:
python $TESS_CODE_DIR/upstream/tica_prefetch.py --sector 102 --cam 3 --ccd 3 --orbit o1a

Every session:
cd /home/$USER/work/TESS/photometry      (where the code lives)
source setup_vars_hpcc.sh                (DATA_DIR/RAW_CACHE on /lustre/scratch)

Whole sector, chained jobs per half-orbit (upstream array -> orbit_pipeline):
bash run_sector_hpcc.sh <sector> [ccd-array]     e.g. 102, or 102 10 for cam3 ccd3

Or step by step:
sbatch --array=0-15%4 upstream/tica_upstream.sbatch <sector> <orbit>
cd $WORK_DIR && sbatch $TESS_CODE_DIR/orbit_pipeline.sbatch <sector> <orbit>

Re-plot already processed orbits (light curves exist) and save each orbit:
sbatch -p nocona -c 8 -t 1-00:00:00 --wrap "bash $TESS_CODE_DIR/reanalyze_orbit.sh <sector> o1a o1b o2a o2b"

To keep using the ISIS products instead, source ~/setup_vars.sh as before
(PHOT_BACKEND then defaults to isis).

LOCAL (WSL) RUNS:

In an Ubuntu (WSL) terminal, from the repository:

source setup_vars_local.sh      (replaces ~/setup_vars.sh: conda env "tess",
                                 DATA_DIR/WORK_DIR on D:, local sbatch/squeue)

1) Upstream processing. On HPCC the image-subtraction products under DATA_DIR
   were made by someone else's ISIS processing. upstream/ rebuilds the same
   products from MAST TICA FFIs (download -> subtract -> delete raw):

   bash upstream/tica_upstream.sh <sector> [cam|all] [ccd|all] [orbit|all]
   (cam/ccd/orbit default to all; stops if the D: drive would fill up)

   $DATA_DIR/s00XX/camC-ccdD/ref.fits, ref_epsf.fits
   $DATA_DIR/s00XX/camC-ccdD/oXX/rms_list_tica
   $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/dates
   $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/images.tar  (conv_*, interp_*, bkg_*)

2) The usual orbit pipeline, run from WORK_DIR (the local sbatch runs each job
   to completion; photometry uses phot_slice.py instead of ISIS phot2):

   cd $WORK_DIR
   bash $TESS_CODE_DIR/orbit_pipeline.sbatch <sector> [orbit|all]
   (orbit defaults to all half-orbits found under DATA_DIR)

   After an orbit is analysed its images.tar files are deleted (PURGE_IMAGES=1,
   the local default) and the orbit is marked .purged so reruns skip it.

3) Or both steps for a whole sector, one half-orbit at a time (disk holds only
   one half-orbit of difference images at once):

   bash run_sector_local.sh <sector> [cam|all] [ccd|all]

Lossy-compression test: prefix both steps with COMPRESS_Q=4 to work in
parallel trees tica_q4/ and photometry_q4/ built from fpack-compressed FFIs.
Settings are listed in setup_vars_local.sh.

PRE PIPELINE:

First: run make_dir.sbatch to create directories

sbatch make_dir.sbatch <sector>

Second: create daily rms images using make_rms_sbatch.sh

sbatch make_rms_sbatch.sh <sector>

PIPELINE:

1) run source.sbatch for the given <sector>:

sbatch source.sbatch <sector> (for entire sector)

sbatch --array=0-0 source.sbatch <sector> <cam> <ccd> <orbit> (for specific cam/ccd/orbit)

2) run photometry with submit_all.sbatch

3) run image_untar.sbatch to untar and produce the slice count image products

sbatch images_untar.sbatch 101

4) perform analysis and make plots with analysis.sbatch

sbatch analysis.sbatch 101 
