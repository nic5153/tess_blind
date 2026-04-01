#!/bin/bash
#SBATCH --job-name=copyclean
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%j
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%j

source ~/setup_vars.sh

python /lustre/research/mfausnau/software/tess_image_sub/phot_scripts/concat_lc.py --sector 76 --cam 4 --ccd 2 --lcdir lcGRB --override --delete
