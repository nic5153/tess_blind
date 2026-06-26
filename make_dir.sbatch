#!/usr/bin/env bash
#SBATCH --job-name="make_sector_dirs"
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%j
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-0:10:00

set -eo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <sector>"
    echo "Example: sbatch $0 101"
    exit 1
fi

sector=$1
sector_dir="sector${sector}"

echo "Creating local directory tree: $PWD/$sector_dir"

mkdir -p "$sector_dir"

for cam in 1 2 3 4; do
    for ccd in 1 2 3 4; do
        camccd_dir="$sector_dir/cam${cam}_ccd${ccd}"

        mkdir -p "$camccd_dir"

        for orbit in o1a o1b o2a o2b; do
            mkdir -p "$camccd_dir/$orbit"
        done
    done
done

echo "Done."
find "$sector_dir" -maxdepth 2 -type d | sort
EOF

bash -n make_sector_dirs.sbatch
