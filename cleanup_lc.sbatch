#!/usr/bin/env bash
#SBATCH --job-name="cleanup_lc"
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%A_%a
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%A_%a
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
dtarget="${DATA_DIR}/${sectoruse}/cam${cam}-ccd${ccd}"

echo "cleanup_lc: sector=$sectoruse cam=$cam ccd=$ccd lcdir=$lcdir  task=$SLURM_ARRAY_TASK_ID"

if [[ ! -d "$dtarget" ]]; then
    echo "ERROR: Target directory not found: $dtarget"
    exit 1
fi

# Build flat list of all slice directories across all orbits
mapfile -t all_slices < <(ls -d "$dtarget"/o??/slice* 2>/dev/null)
nelements=${#all_slices[@]}

if [[ $nelements -eq 0 ]]; then
    echo "No slice directories found, nothing to do."
    exit 0
fi

echo "Total slice directories: $nelements"

clean_lc_contents() {
    local dir=$1
    if [[ ! -d "$dir" ]]; then
        return 0
    fi
    if ! find "$dir" -mindepth 1 -maxdepth 1 -exec rm -rf -- {} +; then
        echo "ERROR: could not fully remove contents of $dir" >&2
        return 1
    fi
    if [[ -n "$(find "$dir" -mindepth 1 -maxdepth 1 -print -quit)" ]]; then
        echo "ERROR: contents remain in $dir after cleanup" >&2
        return 1
    fi
}

for (( ii=1; ii<=nelements; ii++ )); do
    if [[ $(( ii % n_tasks )) -ne $SLURM_ARRAY_TASK_ID ]]; then
        continue
    fi

    slice_dir="${all_slices[$((ii-1))]}"

    clean_lc_contents "$slice_dir/$lcdir"
    clean_lc_contents "$slice_dir/bkg_phot/$lcdir"
done

echo "Task $SLURM_ARRAY_TASK_ID done."
