#!/bin/bash
#SBATCH --job-name=untar_slices
#SBATCH --array=0-15
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --output=/lustre/scratch/nimcclur/rms_testing/logs/%x.o%A_%a
#SBATCH --error=/lustre/scratch/nimcclur/rms_testing/logs/%x.e%A_%a

if [ -z "$1" ]; then
    echo "Usage: sbatch untar_all_slices_live.sbatch <sector_number>"
    exit 1
fi

SECTOR_NUM=$1
SECTOR_DIR=$(printf "s%04d" $SECTOR_NUM)

cams=(1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4)
ccds=(1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4)

CAM=${cams[$SLURM_ARRAY_TASK_ID]}
CCD=${ccds[$SLURM_ARRAY_TASK_ID]}

BASE_PATH="/lustre/research/mfausnau/data/tica/${SECTOR_DIR}/cam${CAM}-ccd${CCD}"

find "$BASE_PATH" -name "images.tar" -execdir echo "Extracting in {}" \; -execdir tar -xf images.tar \;

if [[ $? -eq 0 ]]; then
    echo "------------------------------------------------------------"
    echo "SUCCESS: All images extracted for Cam${CAM}-CCD${CCD}."
    echo "------------------------------------------------------------"
else
    echo "ERROR: Extraction failed for one or more slices in Cam${CAM}-CCD${CCD}." >&2
    exit 1
fi
