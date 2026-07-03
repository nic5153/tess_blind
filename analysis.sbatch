#!/usr/bin/env bash

#SBATCH --job-name=cleanANDanalyze
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j 
#SBATCH --partition=nocona 
#SBATCH --nodes 1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1
#SBATCH --time=0-24:00:00

sector=$1
orbit=${2:-}
suse=$(printf "%04d"  $1)
coordfile=coords.txt
indir=lc
RUN_CLEANED=${RUN_CLEANED:-0}
RUN_MATCHBKG=${RUN_MATCHBKG:-1}
ANALYSIS_N=${ANALYSIS_N:-${SLURM_CPUS_PER_TASK:-1}}

if [[ -n "$orbit" && ! "$orbit" =~ ^o[0-9][ab]$ ]]
then
  echo "ERROR: orbit must look like o1a, o1b, o2a, o2b"
  exit 1
fi

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
  if [[ -n "$orbit" ]]
  then
    match_files=()
    for cam in {1..4}
    do
      for ccd in {1..4}
      do
        for rmsfile in cam${cam}_ccd${ccd}/${orbit}/rms_day_*.fits
        do
          [[ -e "$rmsfile" ]] || continue
          rmsbase=$(basename "$rmsfile")
          day=${rmsbase#rms_day_}
          day=${day%.fits}
          for lcfile in cam${cam}_ccd${ccd}/lc/lc_outcatrms_day_${day}_*.txt
          do
            [[ -e "$lcfile" ]] && match_files+=("$lcfile")
          done
        done
      done
    done

    if [[ ${#match_files[@]} -gt 0 ]]
    then
      for ((i=0; i<${#match_files[@]}; i+=500))
      do
        python /home/sarchast/tesstools/matchbkgtolc.py "${match_files[@]:i:500}" || echo "WARNING: matchbkgtolc.py failed for one ${orbit} chunk; continuing with in-plot time matching."
      done
    else
      echo "WARNING: no ${orbit} light curves found for matchbkgtolc.py; continuing."
    fi
  else
    python /home/sarchast/tesstools/matchbkgtolc.py cam*_ccd*/lc/lc_outcatrms_* || echo "WARNING: matchbkgtolc.py failed; continuing with in-plot time matching."
  fi
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

      ## Clear stale per-cam/ccd plots before generating current plots. The
      ## routed single/multiple/flag folders are rebuilt below from this folder,
      ## so leaving old JPGs here can resurrect candidates that are no longer in
      ## the current phot.data.
      rm -rf "cam${cam}_ccd${ccd}/plots"

      plot_args=(--sector "$sector" --photfile "./cam${cam}_ccd${ccd}/phot.data")
      if [[ -n "$orbit" ]]
      then
        plot_args+=(--orbit "$orbit")
      fi
      if [[ "$ANALYSIS_N" -gt 1 ]]
      then
        plot_args+=(--N "$ANALYSIS_N")
      else
        plot_args+=(--singleproc)
      fi
      python ../complicatedplot.py "${plot_args[@]}"
    fi
  done
done

rm -rf single multiple flag
mkdir -p single multiple flag

move_plot_and_metadata()
{
  local plotfile=$1
  local destdir=$2
  local plotbase
  local metadata

  plotbase=$(basename "$plotfile")
  metadata="${plotfile}.metadata.csv"
  mv "$plotfile" "$destdir/"
  if [[ -e "$metadata" ]]
  then
    mv "$metadata" "$destdir/${plotbase}.metadata.csv"
  fi
}

## Legacy cleaned plot routing:
## mv cam*/plots/flag*cleaned*jpg flag/
## mv cam*/plots/single*outcatrms*jpg single/
## mv cam*/plots/outcatrms*jpg multiple/
if [[ -n "$orbit" ]]
then
  for cam in {1..4}
  do
    for ccd in {1..4}
    do
      for rmsfile in cam${cam}_ccd${ccd}/${orbit}/rms_day_*.fits
      do
        [[ -e "$rmsfile" ]] || continue
        rmsbase=$(basename "$rmsfile")
        day=${rmsbase#rms_day_}
        day=${day%.fits}
        for plotfile in cam${cam}_ccd${ccd}/plots/*outcatrms_day_${day}_*.jpg
        do
          [[ -e "$plotfile" ]] || continue
          plotbase=$(basename "$plotfile")
          if [[ "$plotbase" == flag*outcatrms* ]]
          then
            move_plot_and_metadata "$plotfile" flag
          elif [[ "$plotbase" == single*outcatrms* ]]
          then
            move_plot_and_metadata "$plotfile" single
          else
            move_plot_and_metadata "$plotfile" multiple
          fi
        done
      done
    done
  done
else
  for plotfile in cam*/plots/flag*outcatrms*jpg
  do
    [[ -e "$plotfile" ]] || continue
    move_plot_and_metadata "$plotfile" flag
  done
  for plotfile in cam*/plots/single*outcatrms*jpg
  do
    [[ -e "$plotfile" ]] || continue
    move_plot_and_metadata "$plotfile" single
  done
  for plotfile in cam*/plots/*outcatrms*jpg
  do
    [[ -e "$plotfile" ]] || continue
    move_plot_and_metadata "$plotfile" multiple
  done
fi

echo "Plot page generation skipped (no burst table)."
echo "Plots remain in cam*/plots and the single/ multiple/ flag/ folders."
