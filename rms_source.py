import numpy as np 
import pandas as pd
import argparse 
import subprocess
import os
from astroquery.vizier import Vizier
from tess_stars2px import tess_stars2px_reverse_function_entry
from tqdm import tqdm
from astropy.coordinates import SkyCoord
from astropy.io import fits
import astropy.units as u
import warnings
from pathlib import Path

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

def env_float(name, default):
    try:
        return float(os.environ.get(name, default))
    except ValueError:
        return default

def env_int(name, default):
    try:
        return int(os.environ.get(name, default))
    except ValueError:
        return default

SKIP_BAD_RMS = os.environ.get("RMS_SOURCE_SKIP_BAD_RMS", "1") != "0"
MAX_OUTCAT_ROWS = env_int("RMS_SOURCE_MAX_OUTCAT_ROWS", 50000)
MAX_RMS_MEDIAN = env_float("RMS_SOURCE_MAX_MEDIAN", 1.0e5)
MAX_RMS_P99 = env_float("RMS_SOURCE_MAX_P99", 1.0e8)
MAX_RMS_OVER_MEDIAN = env_float("RMS_SOURCE_MAX_OVER_MEDIAN", 1.0e5)
MIN_FINITE_FRACTION = env_float("RMS_SOURCE_MIN_FINITE_FRACTION", 0.5)

def outcat_to_fits(outcat_file):
    outcat_path = Path(outcat_file)
    if outcat_path.name.startswith("outcat"):
        fits_name = outcat_path.name[len("outcat"):].replace(".tsv", ".fits")
    else:
        fits_name = outcat_path.name.replace(".tsv", ".fits")
    return outcat_path.with_name(fits_name)

def count_outcat_rows(outcat_file):
    with open(outcat_file, "rb") as infile:
        return max(sum(1 for _ in infile) - 1, 0)

def bad_rms_reasons(outcat_file):
    if not SKIP_BAD_RMS:
        return []

    outcat_file = Path(outcat_file)
    fits_file = outcat_to_fits(outcat_file)
    reasons = []

    try:
        nrows = count_outcat_rows(outcat_file)
    except OSError as exc:
        return [f"cannot read outcat ({exc})"]

    if MAX_OUTCAT_ROWS > 0 and nrows > MAX_OUTCAT_ROWS:
        reasons.append(f"outcat_rows={nrows}>{MAX_OUTCAT_ROWS}")

    if not fits_file.exists():
        reasons.append(f"missing FITS {fits_file}")
        return reasons

    try:
        with fits.open(fits_file, memmap=True) as hdul:
            data = np.asarray(hdul[0].data, dtype=float)
            finite = np.isfinite(data)
            finite_fraction = float(finite.mean())
            vals = data[finite]
            if vals.size == 0:
                reasons.append("no finite RMS pixels")
                return reasons

            median = float(np.nanmedian(vals))
            p99 = float(np.nanpercentile(vals, 99))
            maxval = float(np.nanmax(vals))
    except Exception as exc:
        reasons.append(f"cannot read FITS ({exc})")
        return reasons

    if finite_fraction < MIN_FINITE_FRACTION:
        reasons.append(f"finite_fraction={finite_fraction:.3f}<{MIN_FINITE_FRACTION}")
    if (not np.isfinite(median)) or median <= 0:
        reasons.append(f"median={median}")
    elif median > MAX_RMS_MEDIAN:
        reasons.append(f"median={median:.6g}>{MAX_RMS_MEDIAN:.6g}")
    if np.isfinite(p99) and p99 > MAX_RMS_P99:
        reasons.append(f"p99={p99:.6g}>{MAX_RMS_P99:.6g}")
    if np.isfinite(maxval) and np.isfinite(median) and median > 0:
        ratio = maxval / median
        if ratio > MAX_RMS_OVER_MEDIAN:
            reasons.append(f"max_over_median={ratio:.6g}>{MAX_RMS_OVER_MEDIAN:.6g}")

    return reasons

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
tsvfiles = []

for rmsf in args.rms:
    skip_reasons = bad_rms_reasons(rmsf)
    if skip_reasons:
        print(f"SKIP bad RMS source catalog {rmsf}: {'; '.join(skip_reasons)}")
        continue

    rmsdata = np.loadtxt(rmsf, delimiter='\t', usecols=(1,2,4,5), skiprows=1, dtype=[('col','f8'),('row','f8'),('r1','f8'),('r2','f8')])
    if rmsdata.ndim == 0:
        rmsdata = np.array([rmsdata], dtype=rmsdata.dtype)
    rmsfiles.append(rmsdata)
    fitsfiles.append(str(outcat_to_fits(rmsf)))
    tsvfiles.append(rmsf)

objcoords = []
fitsrms = []
tsvfiles = np.array(tsvfiles)
for ff in fitsfiles:
    with fits.open(ff) as hdul:
        fitsrms.append(np.transpose(hdul[0].data))

candidate_id = 0
for rms in rmsfiles:
    for i,src in tqdm(enumerate(rms),total=len(rms)):
        this_candidate_id = candidate_id
        candidate_id += 1
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
        finite_rms = np.isfinite(srcmeasinrms)
        if not finite_rms.any():
            continue
        max_rms = np.nanmax(srcmeasinrms)
        if not np.isfinite(max_rms):
            continue
        other_rms = srcmeasinrms[finite_rms & (srcmeasinrms != max_rms)]
        if (other_rms.size > 0) and ((max_rms/other_rms) < 5).any():
            significant = False
        if (ndets==1) and significant:
            best_rms = np.where(srcmeasinrms==max_rms)[0]
            if best_rms.size == 0:
                continue
            thisfile = tsvfiles[best_rms[0]]
            if args.targetrms is not None:
                if thisfile!=args.targetrms:
                   continue
            ra, dec, scinfo  = tess_stars2px_reverse_function_entry(args.sector, args.cam, args.ccd, src['col'], src['row'], scInfo=scinfo)
            objcoords.append((ra,dec))
            outcat.append((ra,dec,src['col'],src['row'],icol,irow,f"{src['col']+1}\t{src['row']+1}\t{src['r1']}\t{src['r2']}\n",Path(thisfile).stem+f"_cam{args.cam}_ccd{args.ccd}_",this_candidate_id))

objcoords = np.array(objcoords,dtype=[('ra','f8'),('dec','f8')])
print(objcoords.shape)
filerows = []
with open("objphot.tsv","w") as f:
    for ra,dec,fcol,frow,icol,irow,ostring,oname,candidate_id in tqdm(outcat):
        filerows.append(f"{fcol:8.4f}\t{frow:8.4f}\t{icol:d}\t{irow:d}\tlc/lc_{oname}cand{candidate_id}.txt\t1\n")
        f.write(ostring)
    
with open("phot.data","w") as f:
    for o in filerows:
        f.write(o)
