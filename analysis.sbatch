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
indir=lc
RUN_CLEANED=${RUN_CLEANED:-0}
RUN_MATCHBKG=${RUN_MATCHBKG:-1}

source /home/${USER}/setup_vars.sh

cd sector${sector}

## Clean phot file
python ../cleanphotfile.py cam*_ccd*/phot.data

## Legacy cleaned-light-curve workflow.
## Keep this available because it defines the original flux calibration and
## cleaned-column science path, but the current analysis run plots raw lc files.
if [[ "$RUN_CLEANED" == "1" ]]
then
  ## Make coords.txt
  python /home/sarchast/tesstools/phot2coords.py cam*_ccd*/phot.data
fi

## Only get times with bkg and source. This is still useful for the raw plots,
## but complicatedplot.py also aligns by time so this step is not fatal.
if [[ "$RUN_MATCHBKG" == "1" ]]
then
  python /home/sarchast/tesstools/matchbkgtolc.py cam*_ccd*/lc/lc_outcatrms_* || echo "WARNING: matchbkgtolc.py failed; continuing with in-plot time matching."
fi

for cam in {1..4}
do
  for ccd in {1..4}
  do
    if [[ -s cam${cam}_ccd${ccd}/phot.data ]] 
    then 
      echo "check: " $sector $suse $cam $ccd $coordfile $indir
      
      ## Legacy cleaned-light-curve workflow. Leave commented science path in
      ## place, but do not require it for the no-cleaned raw lc analysis.
      if [[ "$RUN_CLEANED" == "1" ]]
      then
        cd cam${cam}_ccd${ccd}

        fluxcal.py --photfile ./phot.data --infiles /lustre/research/mfausnau/data/tica/s"$suse"/"cam$cam""-ccd$ccd"/ref.fits 
        mv ref_gaussian_psf_flux.txt "$indir"/ref_gaussian_psf_flux.txt

        python ../../clean_handmade_lc.py  $indir/lc_outcatrms_* --fluxcal "$indir"/ref_gaussian_psf_flux.txt --metafile coords.txt #--decimal

        cd ..
      fi

      python ../complicatedplot.py --sector $sector --photfile ./cam${cam}_ccd${ccd}/phot.data --singleproc
    fi
  done
done

rm -rf single multiple flag
mkdir -p single multiple flag

## Legacy cleaned plot routing:
## mv cam*/plots/flag*cleaned*jpg flag/
## mv cam*/plots/single*outcatrms*jpg single/
## mv cam*/plots/outcatrms*jpg multiple/
mv cam*/plots/flag*outcatrms*jpg flag/ 2>/dev/null || true
mv cam*/plots/single*outcatrms*jpg single/ 2>/dev/null || true
mv cam*/plots/*outcatrms*jpg multiple/ 2>/dev/null || true

echo "Plot page generation skipped (no burst table)."
echo "Plots remain in cam*/plots and the single/ multiple/ flag/ folders."
