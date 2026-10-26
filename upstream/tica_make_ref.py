"""Build the reference image for one sector/cam/ccd (the role ISIS ref.fits
played on HPCC), plus the star list and per-tile PSFs used for PSF matching.

Follows the MIT TESS-transients convention (Fausnaugh et al.): one reference
per sector and CCD, the median of the 20 lowest-background FFIs from the first
~7 days of the sector. Here: REF_NCAND quality-0 candidates spread evenly over
one half-orbit (o1a by default), of which the REF_NFRAMES lowest-background
frames are background-subtracted, registered onto the middle one's grid, and
median combined.

Outputs in $DATA_DIR/s00XX/camC-ccdD/:
    ref.fits         background-subtracted reference (WCS of the grid frame)
    ref_stars.txt    isolated, unsaturated stars used to measure PSFs
    ref_psf.fits     per-tile reference PSFs (NTILE*NTILE x KSIZE x KSIZE)
    ref_frames.txt   frames combined into the reference
"""
import argparse
import os
import shutil

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.spatial import cKDTree

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--orbit", required=True, help="half-orbit to draw frames from, e.g. o1a")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

outdir = tc.ccd_dir(args.sector, args.cam, args.ccd)
ref_path = os.path.join(outdir, "ref.fits")
if os.path.exists(ref_path) and not args.force:
    print(f"Reference exists: {ref_path} (use --force to rebuild)")
    raise SystemExit(0)

nframes = tc.env_int("REF_NFRAMES", 20)
ncand = max(nframes, tc.env_int("REF_NCAND", 60))
ntile = tc.env_int("NTILE", 4)
ksize = tc.env_int("KSIZE", 13) | 1
max_stars = tc.env_int("MAX_STARS_PER_TILE", 40)
iso_radius = tc.env_float("ISO_RADIUS", 6.0)
q = tc.compress_q()
keep_raw = tc.env_str("KEEP_RAW", "0") == "1"

rows = [r for r in tc.read_manifest(args.sector, args.cam, args.ccd) if r["orbit"] == args.orbit]
if not rows:
    raise SystemExit(f"ERROR: no frames for {args.orbit} in manifest")
rows.sort(key=lambda r: r["cadence"])
picks = np.unique(np.linspace(0, len(rows) - 1, min(ncand, len(rows))).astype(int))
candidates = [rows[i]["name"] for i in picks]

cache = os.path.join(tc.raw_cache(), f"s{args.sector:04d}", f"cam{args.cam}-ccd{args.ccd}", "ref")
print(f"Reference: {len(candidates)} candidates from {args.orbit}, keeping {nframes}"
      + (f", compression q={q}" if q else ""))

# Pass 1: quality flag and background level of every candidate.
levels = []
for name in candidates:
    # take=False: a prefetched frame stays in the pool for its slice job too
    path = tc.fetch_frame(args.sector, args.cam, args.ccd, name, cache, take=False)
    data, hdr, _, _ = tc.load_frame(path, q)
    if int(hdr.get("QUAL_BIT", 0)) != 0:
        print(f"  skip {name}: QUAL_BIT={hdr.get('QUAL_BIT')}")
        continue
    # Frames with missing (non-finite) pixels look deceptively low-background
    # and break Background2D; TICA does not always flag them.
    finite = tc.finite_fraction(data)
    if finite < tc.env_float("REF_MIN_FINITE", 0.99):
        print(f"  skip {name}: only {100 * finite:.1f}% of science pixels are finite")
        continue
    sci = data[tc.SCI_Y, tc.SCI_X][::8, ::8]
    _, med, _ = sigma_clipped_stats(sci, sigma=3.0)
    levels.append((med, name))
if len(levels) < 3:
    raise SystemExit("ERROR: fewer than 3 usable reference candidates")
levels.sort()
chosen = sorted(name for _, name in levels[:nframes])

# Pass 2: background-subtract, register onto the middle chosen frame, median combine.
grid_name = chosen[len(chosen) // 2]
_, grid_hdr, _, _ = tc.load_frame(os.path.join(cache, grid_name), q)
grid_wcs = tc.frame_wcs(grid_hdr)
stack = []
for name in chosen:
    try:
        data, hdr, _, _ = tc.load_frame(os.path.join(cache, name), q)
        bkg, _, _ = tc.science_background(data)
        dx, dy = tc.wcs_shift(grid_wcs, tc.frame_wcs(hdr))
        stack.append(tc.register(data - bkg, dx, dy))
    except Exception as exc:  # one bad frame should not sink the reference
        print(f"  skip {name}: {type(exc).__name__}: {exc}")
        continue
    print(f"  {name}: shift=({dx:+.3f},{dy:+.3f}) px")
if len(stack) < 3:
    raise SystemExit(f"ERROR: only {len(stack)} usable reference frames")
ref = np.median(np.array(stack), axis=0).astype(np.float32)
del stack
mask = np.zeros(ref.shape, bool)
mask[tc.SCI_Y, tc.SCI_X] = True
ref[~mask] = 0.0

# Isolated, unsaturated stars for PSF measurement.
_, _, std = sigma_clipped_stats(ref[tc.SCI_Y, tc.SCI_X], sigma=3.0)
finder = DAOStarFinder(fwhm=tc.env_float("REF_FWHM", 1.8),
                       threshold=tc.env_float("STAR_NSIGMA", 30.0) * std, exclude_border=True)
found = finder(ref)
if found is None or len(found) == 0:
    raise SystemExit("ERROR: no stars found in reference")
xcol = "x_centroid" if "x_centroid" in found.colnames else "xcentroid"
ycol = "y_centroid" if "y_centroid" in found.colnames else "ycentroid"
x = np.asarray(found[xcol], float)
y = np.asarray(found[ycol], float)
flux = np.asarray(found["flux"], float)
half = ksize // 2 + 2
inside = ((x > tc.SCI_X.start + half) & (x < tc.SCI_X.stop - half - 1)
          & (y > tc.SCI_Y.start + half) & (y < tc.SCI_Y.stop - half - 1))
sat = tc.saturation_level(grid_hdr)
unsat = np.array([
    ref[int(round(yy)) - 1:int(round(yy)) + 2, int(round(xx)) - 1:int(round(xx)) + 2].max() < sat
    for xx, yy in zip(x, y)
])
tree = cKDTree(np.c_[x, y])
isolated = np.ones(len(x), bool)
for i, nbrs in enumerate(tree.query_ball_point(np.c_[x, y], iso_radius)):
    others = [j for j in nbrs if j != i]
    if others and np.max(flux[others]) > 0.05 * flux[i]:
        isolated[i] = False
keep = inside & unsat & isolated & (flux > 0)
x, y, flux = x[keep], y[keep], flux[keep]

stars, psfs = [], []
for ty, tx, ys, xs in tc.tiles(ntile):
    sel = np.where((x >= xs.start) & (x < xs.stop) & (y >= ys.start) & (y < ys.stop))[0]
    sel = sel[np.argsort(flux[sel])[::-1][:max_stars]]
    psf, _ = tc.stamp_stack(ref, x[sel], y[sel], ksize)
    if psf is None:
        psf = np.zeros((ksize, ksize))
        psf[ksize // 2, ksize // 2] = 1.0
        print(f"  tile ({ty},{tx}): only {len(sel)} stars, using delta PSF")
    psfs.append(psf)
    for i in sel:
        stars.append((x[i], y[i], ty * ntile + tx, flux[i]))
print(f"  {len(stars)} PSF stars across {ntile}x{ntile} tiles")

hdr = grid_hdr.copy()
hdr["REFNFRM"] = (len(chosen), "frames median-combined into reference")
hdr["REFORBIT"] = (args.orbit, "half-orbit the reference frames came from")
hdr["BKGSUB"] = (True, "reference is background subtracted")
if q:
    hdr["COMPQ"] = (q, "fpack-style RICE quantize level applied to inputs")
os.makedirs(outdir, exist_ok=True)
fits.PrimaryHDU(ref, hdr).writeto(ref_path + ".tmp", overwrite=True, output_verify="silentfix")
fits.PrimaryHDU(np.array(psfs, np.float32)).writeto(
    os.path.join(outdir, "ref_psf.fits"), overwrite=True)
with open(os.path.join(outdir, "ref_stars.txt"), "w") as f:
    f.write("# x y tile flux  (0-based full-frame pixels)\n")
    for sx, sy, t, fl in stars:
        f.write(f"{sx:.3f} {sy:.3f} {t:d} {fl:.6g}\n")
with open(os.path.join(outdir, "ref_frames.txt"), "w") as f:
    f.write("\n".join(chosen) + "\n")
os.replace(ref_path + ".tmp", ref_path)
print(f"Wrote {ref_path}")

if not keep_raw:
    shutil.rmtree(cache, ignore_errors=True)
