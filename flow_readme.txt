FLOW FOR RMS SOURCE AND PHOTOMETRY

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
