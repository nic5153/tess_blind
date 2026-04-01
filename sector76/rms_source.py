import numpy as np 
import pandas as pd
import argparse 
import subprocess
from astroquery.vizier import Vizier
from tess_stars2px import tess_stars2px_reverse_function_entry
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import warnings

def angsep(ra1deg,dec1deg,ra2deg,dec2deg):
    ra1 = ra1deg*np.pi/180
    dec1 = dec1deg*np.pi/180
    ra2 = ra2deg*np.pi/180
    dec2 = dec2deg*np.pi/180
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sep =  np.arccos(np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1 - ra2))
    except Exception as e:
        if round((np.sin(dec1)*np.sin(dec2) + np.cos(dec1)*np.cos(dec2)*np.cos(ra1 - ra2)),1)==1.0:
            sep = 0
    return 180*sep/np.pi

scinfo = None
parser = argparse.ArgumentParser()
parser.add_argument("--sector",type=int,required=True)
parser.add_argument("--cam",type=int, required=True)
parser.add_argument("--ccd",type=int,required=True)
parser.add_argument("--rms",nargs="+",required=True)
parser.add_argument("--targetrms",type=str)
args = parser.parse_args()

outcat= []
minsep = 2 # pixels
rmsfiles = []
fitsfiles = []

for rmsf in args.rms:
    rmsfiles.append(np.loadtxt(rmsf, delimiter='\t', usecols=(1,2,4,5), skiprows=1, dtype=[('col','f8'),('row','f8'),('r1','f8'),('r2','f8')]))
    fitsfiles.append(rmsf.replace("outcat","").replace(".tsv",".fits"))

objcoords = []
fitsrms = []
tsvfiles = np.array(args.rms)
for ff in fitsfiles:
    with fits.open(ff) as hdul:
        fitsrms.append(np.transpose(hdul[0].data))

for rms in rmsfiles:
    for i,src in tqdm(enumerate(rms),total=len(rms)):
        icol = int(round(src['col']))
        irow = int(round(src['row']))
        ndets = 1
        significant = True
        srcmeasinrms = np.zeros(len(tsvfiles))
        for rmsno,(orms,curfits) in enumerate(zip(rmsfiles,fitsrms)):
            srcmeasinrms[rmsno] = curfits[icol,irow]
            if np.array_equal(orms,rms):
                continue
            if ((((src['col'] - orms['col'])**2 + (src['row'] - orms['row'])**2) < minsep**2).any()):
                ndets+=1
        if ((srcmeasinrms.max()/srcmeasinrms[srcmeasinrms!=srcmeasinrms.max()]) < 5).any():
            significant = False
        if (ndets==1) and significant:
            thisfile = tsvfiles[srcmeasinrms==srcmeasinrms.max()][0]
            if args.targetrms is not None:
                if thisfile!=args.targetrms:
                   continue
            ra, dec, scinfo  = tess_stars2px_reverse_function_entry(args.sector, args.cam, args.ccd, src['col'], src['row'], scInfo=scinfo)
            objcoords.append((ra,dec))
            outcat.append((ra,dec,src['col'],src['row'],icol,irow,f"{src['col']+1}\t{src['row']+1}\t{src['r1']}\t{src['r2']}\n",thisfile.replace(".tsv","")+f"_cam{args.cam}_ccd{args.ccd}_"))

objcoords = np.array(objcoords,dtype=[('ra','f8'),('dec','f8')])
print(objcoords.shape)
uno = 0
filerows = []
with open("objphot.tsv","w") as f:
    for ra,dec,fcol,frow,icol,irow,ostring,oname in tqdm(outcat):
        filerows.append(f"{fcol:8.4f}\t{frow:8.4f}\t{icol:d}\t{irow:d}\tlc/lc_{oname}cand{uno}.txt\t1\n")
        uno+=1
        f.write(ostring)
    
with open("phot.data","w") as f:
    for o in filerows:
        f.write(o)
