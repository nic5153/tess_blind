#!/bin/bash 
#SBATCH --job-name=concatlogs  
#SBATCH --nodes=1             
#SBATCH --ntasks-per-node=1   
#SBATCH --cpus-per-task=1     
#SBATCH --partition=nocona    
# SBATCH --account=mfausnau_grp
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%j
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%j
 
source ~/setup_vars.sh
##### cam1-ccd1
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929526* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd1.o23929526 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929526* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd1.e23929526 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929527* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd1.o23929527 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929527* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd1.e23929527 ; rm $logfile ; done
##### cam1-ccd2
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929529* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd2.o23929529 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929529* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd2.e23929529 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929530* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd2.o23929530 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929530* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd2.e23929530 ; rm $logfile ; done
##### cam1-ccd3
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929531* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd3.o23929531 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929531* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd3.e23929531 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929532* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd3.o23929532 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929532* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd3.e23929532 ; rm $logfile ; done
##### cam1-ccd4
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929679* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd4.o23929679 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929679* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam1-ccd4.e23929679 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929680* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd4.o23929680 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929680* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam1-ccd4.e23929680 ; rm $logfile ; done
##### cam2-ccd1
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929730* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd1.o23929730 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929730* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd1.e23929730 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929731* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd1.o23929731 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929731* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd1.e23929731 ; rm $logfile ; done
##### cam2-ccd2
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929733* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd2.o23929733 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929733* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd2.e23929733 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929734* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd2.o23929734 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929734* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd2.e23929734 ; rm $logfile ; done
##### cam2-ccd3
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929804* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd3.o23929804 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929804* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd3.e23929804 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929805* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd3.o23929805 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929805* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd3.e23929805 ; rm $logfile ; done
##### cam2-ccd4
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929807* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd4.o23929807 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929807* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam2-ccd4.e23929807 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929808* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd4.o23929808 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929808* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam2-ccd4.e23929808 ; rm $logfile ; done
##### cam3-ccd1
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929809* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd1.o23929809 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929809* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd1.e23929809 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929810* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd1.o23929810 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929810* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd1.e23929810 ; rm $logfile ; done
##### cam3-ccd2
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929854* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd2.o23929854 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929854* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd2.e23929854 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929855* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd2.o23929855 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929855* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd2.e23929855 ; rm $logfile ; done
##### cam3-ccd3
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929857* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd3.o23929857 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929857* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd3.e23929857 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929858* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd3.o23929858 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929858* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd3.e23929858 ; rm $logfile ; done
##### cam3-ccd4
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929859* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd4.o23929859 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929859* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam3-ccd4.e23929859 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929860* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd4.o23929860 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929860* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam3-ccd4.e23929860 ; rm $logfile ; done
##### cam4-ccd1
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929894* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd1.o23929894 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929894* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd1.e23929894 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929895* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd1.o23929895 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929895* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd1.e23929895 ; rm $logfile ; done
##### cam4-ccd2
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929896* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd2.o23929896 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929896* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd2.e23929896 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929897* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd2.o23929897 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929897* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd2.e23929897 ; rm $logfile ; done
##### cam4-ccd3
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929898* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd3.o23929898 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929898* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd3.e23929898 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929899* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd3.o23929899 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929899* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd3.e23929899 ; rm $logfile ; done
##### cam4-ccd4
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929900* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd4.o23929900 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929900* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/phot_s0076_cam4-ccd4.e23929900 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.o23929901* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd4.o23929901 ; rm $logfile ; done
for logfile in /lustre/work/nimcclur/TESS/photometry/logs/*.e23929901* ; do cat $logfile >> /lustre/work/nimcclur/TESS/photometry/logs/copyclean_s0076_cam4-ccd4.e23929901 ; rm $logfile ; done
