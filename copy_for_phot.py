import argparse
import numpy as np 
import glob
import pandas as pd 
import requests
import shutil
import os 
import tarfile
from astropy.time import Time
from tqdm import tqdm
t1 = Time.now()
parser = argparse.ArgumentParser()
parser.add_argument("--sector",required=True,type=int)
parser.add_argument("--cam",required=True,type=int)
parser.add_argument("--ccd",required=True,type=int)
args = parser.parse_args()

sourcepath = f"/lustre/research/mfausnau/data/tica/f{args.sector:04}/"
outSECpath = f"{os.environ['HOME']}/scratch/s{args.sector:04}"



cam = args.cam
ccd = args.ccd

copyCAMCCD = True
if copyCAMCCD:
    if not os.path.exists(f"{outSECpath}/cam{cam}-ccd{ccd}"):
        os.mkdir(f"{outSECpath}/cam{cam}-ccd{ccd}")
    print("Copy Slice Folders")
    for sliceFOLDER in tqdm(glob.glob(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/cam{cam}-ccd{ccd}/o??/slice*/")):
        folderLIST = sliceFOLDER.split("/")
        orbit = folderLIST[-3]
        curSLICE = folderLIST[-2]
        dates = np.loadtxt(f"{sliceFOLDER}/dates", dtype=[("file","<U128"),("tjd","f8")])
        if not os.path.exists(f"{outSECpath}/cam{cam}-ccd{ccd}/{orbit}"):
            os.mkdir(f"{outSECpath}/cam{cam}-ccd{ccd}/{orbit}")
            for psf_path in glob.glob(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/cam{cam}-ccd{ccd}/{orbit}/psf*"):
                outpsf = psf_path.split("/")[-1]
                outpsfpath = f"{outSECpath}/cam{cam}-ccd{ccd}/{orbit}/{outpsf}"
                if not os.path.exists(outpsfpath):
                    shutil.copy(psf_path,outpsfpath)
        if not os.path.exists(f"{outSECpath}/cam{cam}-ccd{ccd}/{orbit}/{curSLICE}"):
            shutil.copytree(sliceFOLDER,f"{outSECpath}/cam{cam}-ccd{ccd}/{orbit}/{curSLICE}",symlinks=False,ignore_dangling_symlinks=True)
print("Untar images.tar")
for tarfiles in tqdm(glob.glob(f"{outSECpath}/cam{cam}-ccd{ccd}/o??/slice*/images.tar")):
    with tarfile.TarFile(tarfiles) as tf:
        tf.extractall(path=tarfiles.replace("images.tar",""))
    os.remove(tarfiles) 
print("edit procfiles 1")
for procfile in tqdm(glob.glob(f"{outSECpath}/cam{cam}-ccd{ccd}/o??/slice*/process_config")):
    procpath = procfile.split("/")
    camfolder = procpath[-4]
    orbitfolder = procpath[-3]
    with open(procfile,"r") as f:
        oldtext = f.read()
        newtext = oldtext.replace(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/{camfolder}/{orbitfolder}/slice",f"{outSECpath}/{camfolder}/{orbitfolder}/slice")
        newtext = newtext.replace(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/{camfolder}/{orbitfolder}/phot.data",f"{outSECpath}/{camfolder}/{orbitfolder}/phot.data")
    with open(procfile,"w") as f:
        f.write(newtext)
    
print("edit procfiles 2")
for procfile in tqdm(glob.glob(f"{outSECpath}/cam{cam}-ccd{ccd}/o??/slice*/bkg_phot/process_config")):
    procpath = procfile.split("/")
    camfolder = procpath[-5]
    orbitfolder = procpath[-4]
    with open(procfile,"r") as f:
        oldtext = f.read()
        newtext = oldtext.replace(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/{camfolder}/{orbitfolder}/slice",f"{outSECpath}/{camfolder}/{orbitfolder}/slice")
        newtext = newtext.replace(f"/lustre/research/mfausnau/data/tica/s{args.sector:04}/{camfolder}/{orbitfolder}/phot.data",f"{outSECpath}/{camfolder}/{orbitfolder}/phot.data")
    with open(procfile,"w") as f:
        f.write(newtext)
    
t2 = Time.now()

print(float((t2-t1).to_value('hr')),'hours elapsed')
