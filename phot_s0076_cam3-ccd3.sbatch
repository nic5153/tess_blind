#!/bin/bash
#SBATCH --job-name=phot2
#SBATCH --array=0-49
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%A_%a
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%A_%a

source ~/setup_vars.sh

folders=($(eval "ls -d ~/scratch/s0076/cam3-ccd3/o??/slice*/"))
for i in {1..1907}
do
  if [[ $((i%50)) -eq $SLURM_ARRAY_TASK_ID ]]
  then
    thisfolder=${folders[$i-1]}
    cd $thisfolder
    echo "----------" 
    echo $thisfolder 
    echo "----------" 1>&2
    echo $thisfolder 1>&2
    echo 'Decompress kernel_data if it exists'1>&2
    [[ -e 'kernel_data.tgz' ]] && { tar -xzf kernel_data.tgz; }
    tarstat=$?
    sleep 5
    [[ -e psf_table ]] || { ln -s ../psf_table ; }
    [[ -e process_config ]] || { ln -s ../process_config . ; }
    for p in $(ls ../psf_file*fits) ; do [[ -e ${p:3:30} ]] || { ln -s $p ; } ; done 
    echo 'tar exit status' ${tarstat}
    echo 'Run phot2.csh'
    ${ISIS_DIR}/phot2.csh
    echo 'tar kernel_data.tgz if it does not exist'
    [[ -e 'kernel_data.tgz' ]] || tar -czf kernel_data.tgz {kt_,kc_}*fits
    sleep 5
    [[ -e 'bkg_phot' ]] || { mkdir bkg_phot ; cd bkg_phot; ${PIPELINE_DIR}/setup/make_bkg_phot_dir_em2; cd .. ; } 
    cd bkg_phot
    ln -s ../*.fits .
    [[ -e 'lcGRB' ]] && { rm -rf lcGRB ; } 
    mkdir lcGRB
    echo 'Decompress kernel_data if it exists'1>&2
    [[ -e 'kernel_data.tgz' ]] && { tar -xzf kernel_data.tgz; }
    tarstat=$?
    sleep 5
    [[ -e psf_table ]] || { ln -s ../psf_table ; }
    [[ -e process_config ]] || { ln -s ../process_config . ; }
    for p in $(ls ../psf_file*fits) ; do [[ -e ${p:3:30} ]] || { ln -s $p ; } ; done 
    echo 'tar exit status' ${tarstat}
    echo 'Run phot2.csh'
    ${ISIS_DIR}/phot2.csh
    [[ -e 'kernel_data.tgz' ]] || tar -czf kernel_data.tgz {kt_,kc_}*fits
    [[ -e 'kernel_data.tgz' ]] && { rm -rf kt_*fits kc_*fits; }
    cd ..
    [[ -e 'kernel_data.tgz' ]] && { rm -rf kt_*fits kc_*fits; }

    echo "----------" 1>&2
    echo "----------" 
  fi
done
