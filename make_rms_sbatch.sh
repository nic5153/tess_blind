#!/usr/bin/env bash

#SBATCH --job-name=rms_nocopy
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G

logdir="${LOG_DIR:-/lustre/work/nimcclur/TESS/photometry/logs}"
mkdir -p "$logdir"

SECTOR=$1
CAM=$2
CCD=$3
ORBIT=${4:-}

if [[ -n "$ORBIT" ]]; then
  python /home/nimcclur/work/TESS/photometry/make_rms.py --sector $SECTOR --cam $CAM --ccd $CCD --orbit "$ORBIT"
else
  python /home/nimcclur/work/TESS/photometry/make_rms.py --sector $SECTOR --cam $CAM --ccd $CCD
fi

shopt -s nullglob
logs=(rms_nocopy.o* rms_nocopy.e*)
if [[ ${#logs[@]} -gt 0 ]]; then
  mv -- "${logs[@]}" "$logdir"/
fi
