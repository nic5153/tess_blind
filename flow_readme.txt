FLOW FOR RMS SOURCE AND PHOTOMETRY

HPCC CODE UPDATE SAFETY:

The HPCC photometry directory also contains run products, logs, plot tarballs,
and recovered review files. Avoid switching branches there unless you know the
working tree is clean and disposable.

Safer update path for HPCC:

bash sync_code_from_codex.sh

This fetches GitHub and restores only code/script files from
origin/codex/raw-lc-daily-rms, leaving sector outputs and plot exports in place.

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
