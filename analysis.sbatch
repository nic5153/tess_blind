#!/usr/bin/env bash

#SBATCH --job-name=cleanANDanalyze
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j 
#SBATCH --partition=nocona 
#SBATCH --nodes 1
#SBATCH --cpus-per-task 1 
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00

sector=$1
suse=$(printf "%04d"  $1)
coordfile=coords.txt
indir=lcGRB

source /home/${USER}/setup_vars.sh


## Clean phot file
python ../cleanphotfile.py c*/phot.data
## Make coords.txt
python /home/sarchast/tesstools/phot2coords.py c*/phot.data
## Only get times with bkg and source
python /home/sarchast/tesstools/matchbkgtolc.py c*/lcGRB/outcatrms_*

for cam in {1..4}
do
  for ccd in {1..4}
  do
    if [[ -e cam${cam}_ccd${ccd}/coords.txt ]] 
    then 
      cd cam${cam}_ccd${ccd}
      echo "check: " $sector $suse $cam $ccd $ccordifle $indir
      
      fluxcal.py --photfile ./phot.data --infiles /lustre/research/mfausnau/data/tica/s"$suse"/"cam$cam""-ccd$ccd"/ref.fits 
      mv ref_gaussian_psf_flux.txt "$indir"/ref_gaussian_psf_flux.txt
      
      python ../../clean_handmade_lc.py  $indir/outcatrms_* --fluxcal "$indir"/ref_gaussian_psf_flux.txt --metafile coords.txt #--decimal
      

      cd .. 
      python ../complicatedplot.py --sector $sector --photfile ./cam${cam}_ccd${ccd}/phot.data --singleproc
    fi
  done
done

[[ -e single ]] ; { rm -rf single ; }
[[ -e multiple ]] ; { rm -rf multiple ; }
[[ -e flag ]] ; { rm -rf flag ; }
mkdir single
mkdir multiple
mkdir flag
mv cam*/plots/flag*cleaned*jpg flag/
mv cam*/plots/single*outcatrms*jpg single/
mv cam*/plots/outcatrms*jpg multiple/

echo "Plot page generation skipped (no burst table)."
echo "Plots remain in cam*/plots and the single/ multiple/ flag/ folders."
