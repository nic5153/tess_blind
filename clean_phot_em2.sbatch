#!/usr/bin/env bash
#SBATCH --job-name="clean_phot"
#SBATCH --output=/lustre/research/mfausnau/logs/%x.o%A_%a
#SBATCH --error=/lustre/research/mfausnau/logs/%x.e%A_%a
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-200%50
#SBATCH --time=0-1:30:00

set -eo pipefail

# Number of array tasks — used as modulo for file distribution.
# This array is 1-based (--array=1-N), so n_tasks = SLURM_ARRAY_TASK_MAX directly.
n_tasks=$SLURM_ARRAY_TASK_MAX

duse=$DATA_DIR
dhome=$(pwd)
sector=$1
sectoruse=$(printf "s%04d" "$1")
cam=$2
ccd=$3
lcdir=$4

echo "sector use: $sectoruse  cam: $cam  ccd: $ccd  lcdir: $lcdir"

dtarget="$duse/$sectoruse/cam${cam}-ccd${ccd}"
outbase="$dhome/sector${sector}/cam${cam}_ccd${ccd}"

# Build file list excluding .png and already-cleaned files
mapfile -t infiles < <(
    ls "$outbase/$lcdir/lc_"* 2>/dev/null \
    | grep -v "\.png$" \
    | grep -v "\.cleaned$" \
    | awk -F"/" '{print $NF}'
)
nelements=${#infiles[@]}
echo "Total lc files to process: $nelements"

if [[ $nelements -eq 0 ]]; then
    echo "ERROR: No lc_* files found in $outbase/$lcdir/"
    exit 1
fi

for (( ii=$SLURM_ARRAY_TASK_ID; ii<=$nelements; ii+=$n_tasks )); do
    f="${infiles[$((ii-1))]}"
    echo "Processing: $f  (ii=$ii  task=$SLURM_ARRAY_TASK_ID)"

    lc_main="$outbase/$lcdir/$f"
    lc_bkg="$outbase/bkg_phot/$lcdir/$f"

    # Sort and deduplicate by first column (timestamp), in-place
    sort -nuk1 "$lc_main" > "${lc_main}.tmp" && mv "${lc_main}.tmp" "$lc_main"
    sort -nuk1 "$lc_bkg"  > "${lc_bkg}.tmp"  && mv "${lc_bkg}.tmp"  "$lc_bkg"

    ncheck=$( wc -l < "$lc_main" )
    ncheck2=$( wc -l < "$lc_bkg" )

    printf "%20s %20s\n"   "N_lc_main"  "N_lc_bkg"
    printf "%20s %20s\n"   "$ncheck"    "$ncheck2"

    # Count total lines across all orbit/slice directories in DATA_DIR
    nlines=0
    for o_dir in "$dtarget"/o??; do
        [[ -d "$o_dir" ]] || continue
        for slice_dir in "$o_dir"/slice*; do
            [[ -d "$slice_dir" ]] || continue
            lc_slice="$slice_dir/$lcdir/$f"
            if [[ -e "$lc_slice" ]]; then
                jj=$( wc -l < "$lc_slice" )
                (( nlines += jj )) || true
            fi
        done
    done

    printf "%20s %20s %20s\n" "N_lc_main" "N_lc_bkg" "N_lines_in_slices"
    printf "%20s %20s %20s\n" "$ncheck"   "$ncheck2"  "$nlines"

    if [[ "$ncheck" -ne "$ncheck2" ]]; then
        echo "WARNING: $f — main ($ncheck) and bkg ($ncheck2) line counts differ. Skipping deletion."
        continue
    fi

    if [[ "$ncheck" -ne "$nlines" ]]; then
        echo "WARNING: $f — aggregated ($ncheck) != sum of slices ($nlines). Skipping deletion."
        continue
    fi

    echo "$f passed validation. Removing from slice directories."
    for o_dir in "$dtarget"/o??; do
        [[ -d "$o_dir" ]] || continue
        for slice_dir in "$o_dir"/slice*; do
            [[ -d "$slice_dir" ]] || continue
            rm -f "$slice_dir/$lcdir/$f"
            [[ -d "$slice_dir/bkg_phot" ]] && rm -f "$slice_dir/bkg_phot/$lcdir/$f"
        done
    done

done

# Remove lcdir directories that are now empty.  rmdir is a no-op on non-empty
# dirs, so concurrent array tasks can all attempt this safely.
for o_dir in "$dtarget"/o??; do
    [[ -d "$o_dir" ]] || continue
    for slice_dir in "$o_dir"/slice*; do
        [[ -d "$slice_dir" ]] || continue
        rmdir "$slice_dir/$lcdir"             2>/dev/null || true
        rmdir "$slice_dir/bkg_phot/$lcdir"    2>/dev/null || true
    done
done

echo "Task $SLURM_ARRAY_TASK_ID complete."
