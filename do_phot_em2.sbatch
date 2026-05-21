#!/usr/bin/env bash
#SBATCH --job-name="do_phot"
#SBATCH --output=/lustre/research/mfausnau/logs/%x.o%A_%a
#SBATCH --error=/lustre/research/mfausnau/logs/%x.e%A_%a
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --array=0-99%50
#SBATCH --time=0-0:30:00

set -eo pipefail

# Number of array tasks — used as modulo to distribute slices
n_tasks=$(( SLURM_ARRAY_TASK_MAX + 1 ))

duse=$DATA_DIR
dhome=$(pwd)
sector=$1
sectoruse=$(printf "s%04d" "$1")
cam=$2
ccd=$3
o=$4
lcdir=$5

echo "sector use: $sectoruse  cam: $cam  ccd: $ccd  orbit: $o  lcdir: $lcdir"

mapfile -t folders < <(ls -d "$duse/$sectoruse/cam${cam}-ccd${ccd}/$o/slice"*)
nelements=${#folders[@]}

if [[ $nelements -eq 0 ]]; then
    echo "ERROR: No slice directories found under $duse/$sectoruse/cam${cam}-ccd${ccd}/$o/"
    exit 1
fi

echo "Total slices: $nelements"

for (( ii=1; ii<=nelements; ii++ )); do
    if [[ $(( ii % n_tasks )) -ne $SLURM_ARRAY_TASK_ID ]]; then
        continue
    fi

    slice_dir="${folders[$((ii-1))]}"
    cd "$slice_dir"
    echo "Working in: $(pwd)"

    # Symlink PSF files and table from parent
    for p in ../psf_file*fits; do
        [[ -e "${p:3}" ]] || ln -s "$p"
    done
    [[ -e psf_table ]] || ln -s ../psf_table

    [[ -d "$lcdir" ]] || mkdir -p "$lcdir"

    # Unpack archives if present
    [[ -e images.tar      ]] && tar -xf  images.tar
    [[ -e kernel_data.tgz ]] && tar -xzf kernel_data.tgz

    # Set up bkg_phot dir
    if [[ ! -d bkg_phot ]]; then
        mkdir bkg_phot
        cd bkg_phot
        bash "$PIPELINE_DIR/setup/make_bkg_phot_dir_em2"
        cd ..
    fi

    # Run photometry
    "${ISIS_DIR}/phot2.csh" || true

    # Check for output — phot2.csh may exit 0 even after a segfault
    if ! compgen -G "$lcdir/lc_*" > /dev/null; then
        echo "FAILED: no lc_* output in $slice_dir/$lcdir (phot2.csh likely segfaulted)" >&2
        cd "$dhome"
        continue
    fi

    # Repack images and kernel data; clean up loose fits
    if [[ ! -e images.tar      ]] && ls {interp,conv,bkg}*fits &>/dev/null; then tar -cf  images.tar      {interp,conv,bkg}*fits; fi
    if [[ ! -e kernel_data.tgz ]] && ls {kt_,kc_}*fits         &>/dev/null; then tar -czf kernel_data.tgz {kt_,kc_}*fits;         fi
    rm -f kt_*fits kc_*fits

    # bkg_phot photometry
    cd bkg_phot
    [[ -d "$lcdir" ]] || mkdir -p "$lcdir"
    [[ -e kernel_data.tgz ]] && tar -xzf kernel_data.tgz
    "${ISIS_DIR}/phot2.csh" || true

    # Check for output
    if ! compgen -G "$lcdir/lc_*" > /dev/null; then
        echo "FAILED: no lc_* output in $slice_dir/bkg_phot/$lcdir (phot2.csh likely segfaulted)" >&2
        cd "$dhome"
        continue
    fi
    if [[ ! -e kernel_data.tgz ]] && ls {kt_,kc_}*fits &>/dev/null; then tar -czf kernel_data.tgz {kt_,kc_}*fits; fi
    rm -f kt_*fits kc_*fits
    cd ..

    # Final cleanup of unpacked images
    rm -f {interp,conv,bkg}*fits

    cd "$dhome"
done
