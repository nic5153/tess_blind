#!/usr/bin/python3

import os 
import shutil
import glob
import subprocess
import argparse
import traceback 

parser = argparse.ArgumentParser()
parser.add_argument("sector",type=int)
parser.add_argument("cam",type=int)
parser.add_argument("ccd", type=int)
parser.add_argument("--orbit",nargs="+")
args = parser.parse_args()

username = os.getenv('USER')
duse = f"/home/{username}/scratch/"
dhome = os.getcwd()

sector=args.sector
cam = args.cam
ccd = args.ccd
if args.orbit is None:
    orbits = [o.split('/')[-1] for o in glob.glob(f"{duse}/s{sector:04}/cam{cam}-ccd{ccd}/o??")]
else:
    orbits = args.orbit
print(orbits)
for orbit in orbits:
   print("Copy",f"sector{sector:02}/cam{cam}_ccd{ccd}/phot.data",f"{duse}/s{sector:04}/cam{cam}-ccd{ccd}/{orbit}/phot.data")
   shutil.copyfile(f"sector{sector:02}/cam{cam}_ccd{ccd}/phot.data",f"{duse}/s{sector:04}/cam{cam}-ccd{ccd}/{orbit}/phot.data")
def sbatchARRAYbones(arraysize, cmdlist, jobname):
    contents = f"""#!/bin/bash
#SBATCH --job-name={jobname}
#SBATCH --array=0-49
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --output={os.environ['LOG_DIR']}/%x.o%A_%a
#SBATCH --error={os.environ['LOG_DIR']}/%x.e%A_%a

source ~/setup_vars.sh

folders=($(eval "ls -d ~/scratch/s{sector:04}/cam{cam}-ccd{ccd}/o??/slice*/"))
for i in {{1..{arraysize}}}
do
  if [[ $((i%50)) -eq $SLURM_ARRAY_TASK_ID ]]
  then
    thisfolder=${{folders[$i-1]}}
    cd $thisfolder
    echo "----------" 
    echo $thisfolder 
    echo "----------" 1>&2
    echo $thisfolder 1>&2
"""
    for c in cmdlist:
        contents += "    "+c+"\n"
    contents += """
    echo "----------" 1>&2
    echo "----------" 
  fi
done
"""
    return contents

def sbatchBONES(cmdlist, jobname):
    contents = f"""#!/bin/bash
#SBATCH --job-name={jobname}
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --partition=nocona
#SBATCH --output={os.environ['LOG_DIR']}/%x.o%j
#SBATCH --error={os.environ['LOG_DIR']}/%x.e%j

source ~/setup_vars.sh

"""
    for c in cmdlist:
        contents += c+"\n"
    return contents


arraysize = 0
cmdlist = []
for orbit in orbits:
    os.chdir(f"{duse}/s{sector:04}/cam{cam}-ccd{ccd}/{orbit}/")
    odir = os.getcwd()
    for sliceFOLDER in glob.glob("slice*/"):
        os.chdir(odir+"/"+sliceFOLDER)
        arraysize += 1
        if os.path.exists("lcGRB"):
            shutil.rmtree("lcGRB")
        try:
            os.mkdir("lcGRB")
        except Exception as e:
            print(e)
            print(os.getcwd())
            
        if not os.path.exists("psf_table"):
            if os.path.exists("../psf_table"):
                os.symlink("../psf_table","psf_table")
            else:
                raise Exception("Missing psf table")
            # print("ln -s ../psf_table","psf_table")
        for p in glob.glob("../psf_file*fits"):
            if not os.path.exists(p.split('/')[-1]):
                os.symlink(p,p.split('/')[-1])
                # print("ln -s",p,p.split('/')[-1])
cmdlist = ["""echo 'Decompress kernel_data if it exists'1>&2""","""[[ -e 'kernel_data.tgz' ]] && { tar -xzf kernel_data.tgz; }""","""tarstat=$?""","""sleep 5""","""[[ -e psf_table ]] || { ln -s ../psf_table ; }""","""[[ -e process_config ]] || { ln -s ../process_config . ; }""","""for p in $(ls ../psf_file*fits) ; do [[ -e ${p:3:30} ]] || { ln -s $p ; } ; done ""","""echo 'tar exit status' ${tarstat}""","""echo 'Run phot2.csh'""","""${ISIS_DIR}/phot2.csh""","""echo 'tar kernel_data.tgz if it does not exist'""","""[[ -e 'kernel_data.tgz' ]] || tar -czf kernel_data.tgz {kt_,kc_}*fits""","""sleep 5""","""[[ -e 'bkg_phot' ]] || { mkdir bkg_phot ; cd bkg_phot; ${PIPELINE_DIR}/setup/make_bkg_phot_dir_em2; cd .. ; } ""","""cd bkg_phot""","""ln -s ../*.fits .""","""[[ -e 'lcGRB' ]] && { rm -rf lcGRB ; } ""","""mkdir lcGRB""","""echo 'Decompress kernel_data if it exists'1>&2""","""[[ -e 'kernel_data.tgz' ]] && { tar -xzf kernel_data.tgz; }""","""tarstat=$?""","""sleep 5""","""[[ -e psf_table ]] || { ln -s ../psf_table ; }""","""[[ -e process_config ]] || { ln -s ../process_config . ; }""","""for p in $(ls ../psf_file*fits) ; do [[ -e ${p:3:30} ]] || { ln -s $p ; } ; done ""","""echo 'tar exit status' ${tarstat}""","""echo 'Run phot2.csh'""","""${ISIS_DIR}/phot2.csh""","""[[ -e 'kernel_data.tgz' ]] || tar -czf kernel_data.tgz {kt_,kc_}*fits""","""[[ -e 'kernel_data.tgz' ]] && { rm -rf kt_*fits kc_*fits; }""","""cd ..""","""[[ -e 'kernel_data.tgz' ]] && { rm -rf kt_*fits kc_*fits; }"""]

phot2BATCHfile = sbatchARRAYbones(arraysize, cmdlist, "phot2")
os.chdir(dhome)
with open(f"phot_s{sector:04}_cam{cam}-ccd{ccd}.sbatch","w") as f:
    f.write(phot2BATCHfile)

cmdlist = [f"""python /lustre/research/mfausnau/software/tess_image_sub/phot_scripts/concat_lc.py --sector {sector} --cam {cam} --ccd {ccd} --lcdir lcGRB --override --delete"""]
copycleanBATCHfile = sbatchBONES(cmdlist,"copyclean")
with open(f"copyclean_s{sector:04}_cam{cam}-ccd{ccd}.sbatch","w") as f:
    f.write(copycleanBATCHfile)                    
                    
                
        
