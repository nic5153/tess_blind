"""Build per-tile effective PSFs (photutils EPSFBuilder) from the reference image
and its PSF-star list, for the forced PSF photometry in phot_slice.py (the
psf_table / psf_file*.fits role on HPCC).

Writes $DATA_DIR/s00XX/camC-ccdD/ref_epsf.fits: an (NTILE*NTILE, ny, nx) cube of
oversampled ePSFs (header: OVERSAMP, NTILE). Tiles with fewer than
EPSF_MIN_STARS stars fall back to the ePSF built from all stars.
"""
import argparse
import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table
from photutils.psf import EPSFBuilder, extract_stars

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

ccd_dir = tc.ccd_dir(args.sector, args.cam, args.ccd)
out_path = os.path.join(ccd_dir, "ref_epsf.fits")
if os.path.exists(out_path) and not args.force:
    print(f"ePSF exists: {out_path} (use --force to rebuild)")
    raise SystemExit(0)

ntile = tc.env_int("NTILE", 4)
oversamp = tc.env_int("EPSF_OVERSAMP", 2)
size = tc.env_int("EPSF_SIZE", 15) | 1
min_stars = tc.env_int("EPSF_MIN_STARS", 10)

ref = fits.getdata(os.path.join(ccd_dir, "ref.fits")).astype(float)
stars = np.loadtxt(os.path.join(ccd_dir, "ref_stars.txt"), ndmin=2)
nddata = NDData(ref)


def build(sel):
    cutouts = extract_stars(nddata, Table({"x": stars[sel, 0], "y": stars[sel, 1]}), size=size)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        epsf, _ = EPSFBuilder(oversampling=oversamp, maxiters=10, progress_bar=False)(cutouts)
    return np.asarray(epsf.data, dtype=np.float32)


global_epsf = build(np.ones(len(stars), bool))
cube = []
for t in range(ntile * ntile):
    sel = stars[:, 2].astype(int) == t
    if sel.sum() >= min_stars:
        try:
            cube.append(build(sel))
            continue
        except Exception as exc:  # EPSFBuilder can fail on very sparse tiles
            print(f"  tile {t}: EPSFBuilder failed ({exc}); using global ePSF")
    else:
        print(f"  tile {t}: {sel.sum()} stars; using global ePSF")
    cube.append(global_epsf)

hdr = fits.Header()
hdr["OVERSAMP"] = (oversamp, "ePSF oversampling factor")
hdr["NTILE"] = (ntile, "tiles per side; cube index = ty * NTILE + tx")
hdr["NSTARS"] = (len(stars), "PSF stars from ref_stars.txt")
fits.PrimaryHDU(np.array(cube), hdr).writeto(out_path, overwrite=True)
print(f"Wrote {out_path}: {len(cube)} ePSFs, shape {cube[0].shape}, oversampling {oversamp}")
