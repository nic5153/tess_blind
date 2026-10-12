"""Python stand-in for ISIS phot2.csh: forced PSF photometry (photutils
PSFPhotometry with the per-tile ePSFs from upstream/tica_make_epsf.py) on the
difference images of one slice, at the positions in the orbit's phot.data.

Called by do_phot_em2.sbatch when PHOT_BACKEND=python. Writes, like phot2,

    <slice>/<lcdir>/lc_*.txt            difference-image light curves
    <slice>/bkg_phot/<lcdir>/lc_*.txt   background light curves

one row per frame, in the column layout the downstream scripts read:

    0  time [TJD, exposure start, from the slice's dates file]
    1  difference flux [e-], ISIS sign (reference minus image: brightening < 0)
    2  flux uncertainty [e-]
    3  x, 4  y  (0-based pixel centers on the reference grid)
    5  fit quality (photutils qfit)
    6  sky level at the source [e-] (frame sky, from bkg_<name>)

In bkg_phot the flux columns measure the residual background around the source
in the same difference image: the sigma-clipped median in an annulus
(PHOT_BKG_ANNULUS, default 4-8 px) expressed as the flux a fixed-position PSF
fit would assign to it (median x sum(P)/sum(P^2)). It shows whether a
"brightening" is just local background or scattered-light residuals.

phot.data positions are DAOStarFinder centroids minus 0.5 (do_photutils_extract.py),
so PHOT_XY_OFFSET (default 0.5) is added back before fitting.
"""
import argparse
import io
import os
import sys
import tarfile
import warnings

import numpy as np
from astropy.io import fits
from astropy.table import Table
from photutils.psf import ImagePSF, PSFPhotometry

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "upstream"))
import tica_common as tc  # noqa: E402

parser = argparse.ArgumentParser()
parser.add_argument("--slice-dir", required=True)
parser.add_argument("--lcdir", default="lc")
parser.add_argument("--photfile", default=None, help="default: <slice>/../phot.data")
args = parser.parse_args()

xy_offset = tc.env_float("PHOT_XY_OFFSET", 0.5)
fit_size = tc.env_int("PHOT_FIT_SHAPE", 5) | 1
aper = tc.env_float("PHOT_APERTURE", 2.0)

slice_dir = os.path.abspath(args.slice_dir)
orbit_dir = os.path.dirname(slice_dir)
ccd_dir = os.path.dirname(orbit_dir)
photfile = args.photfile or os.path.join(orbit_dir, "phot.data")

names, xs, ys = [], [], []
with open(photfile) as f:
    for line in f:
        parts = line.split()
        if len(parts) != 6:
            continue
        xs.append(float(parts[0]) + xy_offset)
        ys.append(float(parts[1]) + xy_offset)
        names.append(os.path.basename(parts[4]))
if not names:
    print(f"No sources in {photfile}; nothing to do")
    raise SystemExit(0)
xs, ys = np.array(xs), np.array(ys)

with fits.open(os.path.join(ccd_dir, "ref_epsf.fits")) as hdul:
    epsf_cube = np.array(hdul[0].data, dtype=float)
    oversamp = int(hdul[0].header["OVERSAMP"])
    ntile = int(hdul[0].header["NTILE"])
ref = fits.getdata(os.path.join(ccd_dir, "ref.fits")).astype(float)
tiles = tc.tile_index(xs, ys, ntile)

frames = []
with open(os.path.join(slice_dir, "dates")) as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            frames.append((parts[0], float(parts[1])))


def formal_errors(error, sel, model):
    """Formal flux uncertainty of a fixed-position PSF fit, 1/sqrt(sum(P^2/sigma^2))
    over the fit box. (The fitter's covariance is rescaled by the residuals,
    which is not a noise estimate for a one-parameter linear fit.)"""
    h = fit_size // 2
    out = np.full(sel.sum(), np.nan)
    for k, (x, y) in enumerate(zip(xs[sel], ys[sel])):
        ix, iy = int(round(x)), int(round(y))
        if iy - h < 0 or ix - h < 0 or iy + h >= error.shape[0] or ix + h >= error.shape[1]:
            continue
        yy, xx = np.mgrid[iy - h:iy + h + 1, ix - h:ix + h + 1]
        psf = model.evaluate(xx, yy, 1.0, x, y)
        sig = error[iy - h:iy + h + 1, ix - h:ix + h + 1]
        out[k] = 1.0 / np.sqrt(np.sum(psf ** 2 / sig ** 2))
    return out


bkg_rin, bkg_rout = (float(v) for v in tc.env_str("PHOT_BKG_ANNULUS", "4,8").split(","))


def annulus_flux(image, sel, model):
    """Residual background around each source, as PSF-fit flux (and error)."""
    h = fit_size // 2
    r = int(np.ceil(bkg_rout))
    flux = np.full(sel.sum(), np.nan)
    err = flux.copy()
    for k, (x, y) in enumerate(zip(xs[sel], ys[sel])):
        ix, iy = int(round(x)), int(round(y))
        if iy - r < 0 or ix - r < 0 or iy + r >= image.shape[0] or ix + r >= image.shape[1]:
            continue
        yy, xx = np.mgrid[iy - r:iy + r + 1, ix - r:ix + r + 1]
        rad = np.hypot(xx - x, yy - y)
        ann = image[iy - r:iy + r + 1, ix - r:ix + r + 1][(rad >= bkg_rin) & (rad <= bkg_rout)]
        ann = ann[np.isfinite(ann)]
        if ann.size < 10:
            continue
        med = np.median(ann)
        sig = 1.4826 * np.median(np.abs(ann - med))
        keep = np.abs(ann - med) < 3 * sig if sig > 0 else np.ones(ann.size, bool)
        med = float(np.median(ann[keep]))
        sig = float(np.std(ann[keep]))
        by, bx = np.mgrid[iy - h:iy + h + 1, ix - h:ix + h + 1]
        psf = model.evaluate(bx, by, 1.0, x, y)
        gain = psf.sum() / np.sum(psf ** 2)
        flux[k] = med * gain
        err[k] = sig / np.sqrt(keep.sum()) * gain
    return flux, err


def photometer(image, error, sel, model):
    phot = PSFPhotometry(model, (fit_size, fit_size), aperture_radius=aper)
    init = Table({"x": xs[sel], "y": ys[sel]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = phot(image, error=error, init_params=init)
    flux = np.asarray(res["flux_fit"], float)
    qfit = np.asarray(res["qfit"], float) if "qfit" in res.colnames else np.full(len(flux), np.nan)
    return flux, formal_errors(error, sel, model), qfit


rows = {n: [] for n in names}
bkg_rows = {n: [] for n in names}
with tarfile.open(os.path.join(slice_dir, "images.tar"), "r") as tar:
    for name, tjd in frames:
        with fits.open(io.BytesIO(tar.extractfile(f"conv_{name}").read())) as ch:
            conv = ch[0].data.astype(float)
            sky_level = float(ch[0].header.get("BKGMED", 0.0))
        try:
            with fits.open(io.BytesIO(tar.extractfile(f"bkg_{name}").read())) as bh:
                sky = tc.mesh_to_image(bh[0].data, int(bh[0].header["SKYBOX"]), conv.shape)
        except KeyError:  # slice written without bkg_ maps: flat sky from the header
            sky = np.full_like(conv, sky_level)
        # Poisson variance of the frame [e-]: sky plus stars. The reference is
        # background subtracted (negative between stars in crowded fields), so
        # only its positive part adds star light.
        error = np.sqrt(np.clip(sky + np.clip(ref, 0, None), 1.0, None))

        flux = np.full(len(names), np.nan)
        ferr, qfit = flux.copy(), flux.copy()
        bflux, berr, bqfit = flux.copy(), flux.copy(), flux.copy()
        for t in np.unique(tiles):
            sel = tiles == t
            model = ImagePSF(epsf_cube[t], oversampling=oversamp)
            model.x_0.fixed = True
            model.y_0.fixed = True
            flux[sel], ferr[sel], qfit[sel] = photometer(conv, error, sel, model)
            bflux[sel], berr[sel] = annulus_flux(conv, sel, model)

        ix = np.clip(np.round(xs).astype(int), 0, conv.shape[1] - 1)
        iy = np.clip(np.round(ys).astype(int), 0, conv.shape[0] - 1)
        skyval = sky[iy, ix]
        for i, n in enumerate(names):
            rows[n].append(f"{tjd:.8f} {flux[i]:.6e} {ferr[i]:.6e} {xs[i]:.4f} {ys[i]:.4f} "
                           f"{qfit[i]:.6e} {skyval[i]:.6e}\n")
            bkg_rows[n].append(f"{tjd:.8f} {bflux[i]:.6e} {berr[i]:.6e} {xs[i]:.4f} {ys[i]:.4f} "
                               f"{bqfit[i]:.6e} {skyval[i]:.6e}\n")

for subdir, table in ((args.lcdir, rows), (os.path.join("bkg_phot", args.lcdir), bkg_rows)):
    outdir = os.path.join(slice_dir, subdir)
    os.makedirs(outdir, exist_ok=True)
    for n, lines in table.items():
        with open(os.path.join(outdir, n), "w") as f:
            f.writelines(lines)
print(f"phot_slice: {len(names)} sources x {len(frames)} frames -> {slice_dir}/{args.lcdir}")
