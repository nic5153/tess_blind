#!/bin/bash

source ~/setup_vars.sh
sector=$1
slurmids=()
echo "#!/bin/bash " > cleanlogs.sbatch
echo "#SBATCH --job-name=concatlogs  " >>cleanlogs.sbatch
echo "#SBATCH --nodes=1             ">>cleanlogs.sbatch
echo "#SBATCH --ntasks-per-node=1   ">>cleanlogs.sbatch
echo "#SBATCH --cpus-per-task=1     ">>cleanlogs.sbatch
echo "#SBATCH --partition=nocona    ">>cleanlogs.sbatch
echo "# SBATCH --account=mfausnau_grp">>cleanlogs.sbatch
echo "#SBATCH --output=${LOG_DIR}/%x.o%j">>cleanlogs.sbatch
echo "#SBATCH --error=${LOG_DIR}/%x.e%j">>cleanlogs.sbatch
echo " ">>cleanlogs.sbatch
echo "source ~/setup_vars.sh">>cleanlogs.sbatch

printf -v secstring "%04d" $sector
[[ -e killjobs.sh ]] &&  rm -f killjobs.sh 
# [[ -e prerun.txt ]] && rm -f prerun.txt
# [[ -e postrun.txt ]] && rm -f postrun.txt
# for s in /lustre/research/mfausnau/data/tica/s${secstring}/cam*/o*/slice*/kernel_data.tgz; do du -b ${s} >> prerun.txt ; echo $s; done
workdir=$(pwd)
# The blind search goes over every cam and ccd combo
# for cam in {1..4}
#   do
#   for ccd in {1..4}
#   do 
#     for slice in ~/scratch/s${secstring}/cam${cam}-ccd${ccd}/o??/slice*/
#     do 
#       for kdata in ${slice}/kernel_data.tgz
#       do 
#         out=`wc -l ${slice}/dates`
#         set -- $out
#         nlines=$1
#         minsize=$((nlines*10000))
#         out=`du -b ${kdata}`
#         set -- $out 
#         if [ $1 -lt $minsize ] 
#         then 
#           echo "${kdata} is corrupted"
#           exit 1
#         fi  
#       done 
#     done
#   done
# done
for cam in {1..4}
  do
  for ccd in {1..4}
  do 
    if [[ -e sector${sector}/cam${cam}_ccd${ccd}/phot.data ]] 
    then
      rm -f ~/scratch/s${secstring}/cam${cam}-ccd${ccd}/o??/phot.data 
      cd $workdir 
      python3 do_phot_em2.py ${sector} ${cam} ${ccd} 
      ID=$(sbatch --parsable  phot_s${secstring}_cam${cam}-ccd${ccd}.sbatch)
      slurmids+=($ID)
      echo "##### cam${cam}-ccd${ccd}" >>cleanlogs.sbatch
      echo "for logfile in ${LOG_DIR}/*.o${ID}* ; do cat \$logfile >> ${LOG_DIR}/phot_s${secstring}_cam${cam}-ccd${ccd}.o${ID} ; rm \$logfile ; done" >>cleanlogs.sbatch
      echo "for logfile in ${LOG_DIR}/*.e${ID}* ; do cat \$logfile >> ${LOG_DIR}/phot_s${secstring}_cam${cam}-ccd${ccd}.e${ID} ; rm \$logfile ; done" >>cleanlogs.sbatch 
      echo "scancel ${ID}" >> killjobs.sh 
      ID=$(sbatch --parsable --dependency=afterany:${ID} copyclean_s${secstring}_cam${cam}-ccd${ccd}.sbatch)
      slurmids+=($ID)
      echo "for logfile in ${LOG_DIR}/*.o${ID}* ; do cat \$logfile >> ${LOG_DIR}/copyclean_s${secstring}_cam${cam}-ccd${ccd}.o${ID} ; rm \$logfile ; done" >>cleanlogs.sbatch
      echo "for logfile in ${LOG_DIR}/*.e${ID}* ; do cat \$logfile >> ${LOG_DIR}/copyclean_s${secstring}_cam${cam}-ccd${ccd}.e${ID} ; rm \$logfile ; done" >>cleanlogs.sbatch
      echo "scancel ${ID}" >> killjobs.sh 
    fi
  done
done 
# echo "for s in /lustre/research/mfausnau/data/tica/s${secstring}/cam*/o*/slice*/kernel_data.tgz; do du -b ${s} >> postrun.txt ; done" >> cleanlogs.sbatch
# echo 'diff prerun.txt postrun.txt' >> cleanlogs.sbatch
sbatch --dependency=afterany:$(echo ${slurmids[*]} | tr ' ' :) cleanlogs.sbatch
