#!/bin/bash
#SBATCH --job-name=copy_sector
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=quanah
#SBATCH --array=0-15
#SBATCH --output=logs/%x.o%A_%a
#SBATCH --error=logs/%x.e%A_%a
#SBATCH -t 1-00:00:00


sector=$1
ccd=$((SLURM_ARRAY_TASK_ID%4+1))
j=$((SLURM_ARRAY_TASK_ID/4))
cam=$((j%4+1))


printf -v secstring "%04d" $sector
if [ $cam -eq 1 ] && [ $ccd -eq 1 ] 
  then  
  [[ -e "/home/${USER}/scratch/s${secstring}" ]] || { mkdir "/home/${USER}/scratch/s${secstring}" ; }
else 
  sleep 10
fi

python copy_for_phot.py --sector $sector --cam $cam --ccd $ccd

