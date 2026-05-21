#!/usr/bin/env bash

#SBATCH --job-name=rms_nocopy
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

SECTOR=$1
CAM=$2
CCD=$3

python /home/nimcclur/work/TESS/photometry/rms_daily/make_rms.py --sector $SECTOR --cam $CAM --ccd $CCD

mv rms_nocopy.o* ${LOG_DIR}/
mv rms_nocopy.e* ${LOG_DIR}/
