FLOW FOR RMS SOURCE AND PHOTOMETRY

PRE PIPELINE:

First: run make_dir.sbatch to create directories

sbatch make_dir.sbatch sector76

Second: 

1: create directory for sector

for cam in {1..4} ; do for ccd in {1..4} ; do for orbit o1a o1b o2a o2b ; do mkdir -p cam${cam}_ccd${ccd}/${orbit} ; done ; done ; done

2: copy rms.fits

for cam in {1..4}; do for ccd in {1..4}; do for orbit in o1a o1b o2a o2b; do cp /lustre/research/mfausnau/data/tica/s00**/cam${cam}-ccd${ccd}/${orbit}/rms.fits cam${cam}_ccd${ccd}/${orbit}/; done; done; done

3: inside sector## dir, run do_photutils_extract.sbatch

sbatch do_photutils_extract.sbatch --sector

4: now run rms_source.py inside sector## dir

sbatch rms_source.sbatch --sector

5: run copy_sector.sbatch

sbatch copy_sector.sbatch --sector

6. execute photometry using runscript.sh

./runscript.sh --sector

7. perform analysis, run analysis.sbatch

sbatch analysis.sbatch --sector

