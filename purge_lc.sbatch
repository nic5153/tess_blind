#!/usr/bin/env bash
#SBATCH --job-name="purge_lc"
#SBATCH --output=/lustre/research/mfausnau/logs/%x.o%A_%a
#SBATCH --error=/lustre/research/mfausnau/logs/%x.e%A_%a
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-99%50
#SBATCH --time=0-1:00:00

set -eo pipefail

# Number of array tasks — used as modulo to distribute slices
n_tasks=$(( SLURM_ARRAY_TASK_MAX + 1 ))

sector=$1
cam=$2
ccd=$3
lcdir=$4

sectoruse=$(printf "s%04d" "$sector")
dtarget_base="${DATA_DIR}/${sectoruse}"

echo "purge_lc: sector=$sectoruse cam=$cam ccd=$ccd lcdir=$lcdir  task=$SLURM_ARRAY_TASK_ID"

# Build flat list of all slice directories to process.
# cam/ccd may be "all" (submit_all path) or a specific number.
if [[ "$cam" == "all" ]]; then
    mapfile -t all_slices < <(ls -d "$dtarget_base"/cam*-ccd*/o??/slice* 2>/dev/null)
else
    mapfile -t all_slices < <(ls -d "$dtarget_base/cam${cam}-ccd${ccd}"/o??/slice* 2>/dev/null)
fi

nelements=${#all_slices[@]}

if [[ $nelements -eq 0 ]]; then
    echo "No slice directories found, nothing to do."
    exit 0
fi

echo "Total slice directories: $nelements"

for (( ii=1; ii<=nelements; ii++ )); do
    if [[ $(( ii % n_tasks )) -ne $SLURM_ARRAY_TASK_ID ]]; then
        continue
    fi

    slice_dir="${all_slices[$((ii-1))]}"

    rm -rf "$slice_dir/$lcdir"
    rm -rf "$slice_dir/bkg_phot/$lcdir"
done

echo "Task $SLURM_ARRAY_TASK_ID done."
