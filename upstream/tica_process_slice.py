"""Download and image-subtract one slice, writing the same per-slice products the
HPCC upstream processing produced:

    $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/dates        <name> <STARTTJD>
    $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/images.tar   conv_<name>, interp_<name>, bkg_<name>

dates holds the exposure start time, as on HPCC (clean_handmade_lc.py adds half
an exposure to reach mid-exposure); every image header also carries
STARTTJD/MIDTJD/ENDTJD. bkg_<name> holds two small background maps used by
phot_slice.py: the frame's sky level (primary HDU, SKYBOX-pixel block medians)
and the sky removed from the difference image (DIFFBKG HDU, DIFF_BKG_BOX mesh).

plus bookkeeping used by tica_finalize_orbit.py and the compression tests:

    .quality      <name> <QUAL_BIT> <MIDTJD>
    .stats.tsv    per-frame shift, background, kernel scale, noise, bytes
    .done         marker; the slice is skipped on reruns unless --force

Each frame is resampled onto the reference grid (interp_) and differenced:

  SUB_MODE=reference  conv = ref (x) K - image, K a per-tile PSF-matching kernel
                      KERNEL_METHOD=lstsq      real-space delta-function basis
                                               (Bramich 2008; ISIS-like), default
                      KERNEL_METHOD=photutils  photutils create_matching_kernel
                                               (Fourier; smooths TESS's
                                               undersampled PSF, leaves rings)
  SUB_MODE=adjacent   conv = neighbor - image, the neighboring frame in time
                      (short-timescale transients; no kernel)

The sky is then removed from the difference image with photutils Background2D
(DIFF_BKG_BOX px boxes; QLP-style local background correction), and its rms is
stored as DIFFRMS for noise classification. With COLUMN_CLEAN=1 (default)
column offsets are removed and bleed/smear columns of bright stars are zeroed
(tica_common.clean_columns; header NBADCOL/BADCOLS). conv_ keeps the ISIS sign
(reference minus image) that complicatedplot.py flips.
Raw downloads are deleted once the slice is written unless KEEP_RAW=1.
"""
import argparse
import io
import os
import shutil
import tarfile
import time

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.signal import fftconvolve

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--orbit", required=True)
parser.add_argument("--slice", required=True, help="e.g. slice0003")
parser.add_argument("--force", action="store_true")
args = parser.parse_args()

mode = tc.env_str("SUB_MODE", "reference")
method = tc.env_str("KERNEL_METHOD", "lstsq")
interp_mode = tc.env_str("INTERP_MODE", "header")
ntile = tc.env_int("NTILE", 4)
ksize = tc.env_int("KSIZE", 13) | 1
khalf = tc.env_int("KERNEL_HALF", 2)
diff_box = tc.env_int("DIFF_BKG_BOX", 16)
sky_box = tc.env_int("SKY_BOX", 64)
col_clean = tc.env_str("COLUMN_CLEAN", "1") == "1"
q = tc.compress_q()
keep_raw = tc.env_str("KEEP_RAW", "0") == "1"
if mode not in ("reference", "adjacent"):
    raise SystemExit(f"ERROR: SUB_MODE must be reference or adjacent, not {mode}")
if method not in ("lstsq", "photutils"):
    raise SystemExit(f"ERROR: KERNEL_METHOD must be lstsq or photutils, not {method}")

ccd_dir = tc.ccd_dir(args.sector, args.cam, args.ccd)
slice_dir = os.path.join(ccd_dir, args.orbit, args.slice)
done_path = os.path.join(slice_dir, ".done")
if os.path.exists(done_path) and not args.force:
    print(f"SKIP done {slice_dir}")
    raise SystemExit(0)

orbit_rows = sorted(
    (r for r in tc.read_manifest(args.sector, args.cam, args.ccd) if r["orbit"] == args.orbit),
    key=lambda r: r["cadence"],
)
members = [r for r in orbit_rows if r["slice"] == args.slice]
if not members:
    raise SystemExit(f"ERROR: {args.orbit}/{args.slice} not in manifest")

ref_path = os.path.join(ccd_dir, "ref.fits")
with fits.open(ref_path) as hdul:
    ref = np.array(hdul[0].data, dtype=np.float32)
    ref_hdr = hdul[0].header.copy()
ref_wcs = tc.frame_wcs(ref_hdr)
ref_wcs_hdr = ref_wcs.to_header(relax=True)
sci_mask = np.zeros(ref.shape, bool)
sci_mask[tc.SCI_Y, tc.SCI_X] = True
tile_list = list(tc.tiles(ntile))

cache = os.path.join(tc.raw_cache(), f"s{args.sector:04d}", f"cam{args.cam}-ccd{args.ccd}",
                     args.orbit, args.slice)


def prepare(name):
    path = tc.fetch_frame(args.sector, args.cam, args.ccd, name, cache)
    keep_fz = os.path.join(cache, name + ".fz") if (q is not None and keep_raw) else None
    data, hdr, raw_b, stored_b = tc.load_frame(path, q, keep_fz)
    dx, dy = tc.wcs_shift(ref_wcs, tc.frame_wcs(hdr))
    interp = tc.register(data, dx, dy)
    img = interp.copy()
    img[~sci_mask] = 0.0
    sky_mesh = tc.block_median(interp, sky_box)
    return {
        "name": name, "hdr": hdr, "img": img, "sky_mesh": sky_mesh,
        "interp": interp if interp_mode == "full" else None,
        "dx": dx, "dy": dy, "bkg_med": float(np.median(sky_mesh)),
        "raw_bytes": raw_b, "stored_bytes": stored_b,
        "qual": int(hdr.get("QUAL_BIT", 0)), "tjd": float(hdr["STARTTJD"]),
    }


t0 = time.time()
frames = [prepare(r["name"]) for r in members]

# ---- difference images -------------------------------------------------------
if mode == "reference":
    stars = np.loadtxt(os.path.join(ccd_dir, "ref_stars.txt"), ndmin=2)
    sx, sy, stile = stars[:, 0], stars[:, 1], stars[:, 2].astype(int)
    delta = np.zeros((ksize, ksize))
    delta[ksize // 2, ksize // 2] = 1.0

    if method == "photutils":
        try:
            from photutils.psf_matching import SplitCosineBellWindow, create_matching_kernel
        except ImportError:  # photutils < 3.0
            from photutils.psf.matching import SplitCosineBellWindow, create_matching_kernel
        ref_psfs = fits.getdata(os.path.join(ccd_dir, "ref_psf.fits")).astype(float)
        window = SplitCosineBellWindow(alpha=tc.env_float("WIN_ALPHA", 0.35),
                                       beta=tc.env_float("WIN_BETA", 0.15))
        ref_sums = {t: tc.stamp_stack(ref, sx[stile == t], sy[stile == t], ksize)[1]
                    for t in range(ntile * ntile)}

    def tile_kernel(img, t):
        sel = stile == t
        if method == "lstsq":
            kern = tc.lstsq_kernel(img, ref, sx[sel], sy[sel], half=khalf, ksize=ksize)
            return delta if kern is None or not np.all(np.isfinite(kern)) else kern
        psf, sums = tc.stamp_stack(img, sx[sel], sy[sel], ksize)
        if psf is None:
            return delta
        kern = create_matching_kernel(ref_psfs[t], psf, window=window)
        if not np.all(np.isfinite(kern)) or kern.sum() <= 0:
            return delta
        ok = np.isfinite(sums) & (ref_sums[t] > 0)
        scale = float(np.median(sums[ok] / ref_sums[t][ok])) if ok.sum() >= 3 else 1.0
        return kern / kern.sum() * (scale if 0.5 < scale < 2.0 else 1.0)

    def difference(i):
        img = frames[i]["img"]
        diff = np.zeros_like(ref)
        sums = []
        for t, (_, _, ys, xs) in enumerate(tile_list):
            kern = tile_kernel(img, t)
            sums.append(kern.sum())
            y0, y1 = max(ys.start - ksize, 0), min(ys.stop + ksize, ref.shape[0])
            x0, x1 = max(xs.start - ksize, 0), min(xs.stop + ksize, ref.shape[1])
            model = fftconvolve(ref[y0:y1, x0:x1], kern, mode="same")
            diff[ys, xs] = model[ys.start - y0:ys.stop - y0, xs.start - x0:xs.stop - x0] - img[ys, xs]
        return diff, float(np.median(sums))
else:
    first = orbit_rows.index(members[0])
    if first > 0:
        neighbor_name = orbit_rows[first - 1]["name"]
    elif len(orbit_rows) > 1:
        neighbor_name = orbit_rows[1]["name"]
    else:
        raise SystemExit("ERROR: adjacent mode needs at least two frames in the orbit")
    neighbor0 = (frames[1]["img"] if (first == 0 and len(frames) > 1)
                 else prepare(neighbor_name)["img"])

    def difference(i):
        prev = neighbor0 if i == 0 else frames[i - 1]["img"]
        return prev - frames[i]["img"], 1.0

# ---- write slice products ----------------------------------------------------
os.makedirs(slice_dir, exist_ok=True)
tar_tmp = os.path.join(slice_dir, "images.tar.part")
stats = []


def add_member(tar, name, payload):
    info = tarfile.TarInfo(name)
    info.size = len(payload)
    info.mtime = int(time.time())
    tar.addfile(info, io.BytesIO(payload))


with tarfile.open(tar_tmp, "w") as tar:
    for i, f in enumerate(frames):
        diff, scale = difference(i)
        diff[~sci_mask] = 0.0
        dbkg, _, diff_rms, diff_mesh = tc.science_background(diff, box=diff_box, return_mesh=True)
        diff -= dbkg
        diff[~sci_mask] = 0.0
        badcols = tc.clean_columns(diff) if col_clean else []
        _, _, resid = sigma_clipped_stats(diff[tc.SCI_Y, tc.SCI_X][::4, ::4], sigma=3.0)

        hdr = tc.replace_wcs(f["hdr"], ref_wcs_hdr)
        hdr["SUBMODE"] = (mode, "reference or adjacent-frame subtraction")
        hdr["KMETHOD"] = (method if mode == "reference" else "none", "PSF-matching kernel method")
        hdr["REFFILE"] = ("ref.fits" if mode == "reference" else "neighbor", "subtracted image")
        hdr["SHIFTX"] = (round(f["dx"], 4), "x shift onto reference grid [pix]")
        hdr["SHIFTY"] = (round(f["dy"], 4), "y shift onto reference grid [pix]")
        hdr["BKGMED"] = (round(f["bkg_med"], 3), "median sky level of the frame [e-]")
        hdr["DIFFRMS"] = (round(diff_rms, 3), "median Background2D rms of difference [e-]")
        if q is not None:
            hdr["COMPQ"] = (q, "fpack-style RICE quantize level applied")

        conv_hdr = hdr.copy()
        conv_hdr["DIFFSIGN"] = ("REF-IMG", "ISIS convention: subtracted minus image")
        conv_hdr["KSCALE"] = (round(scale, 5), "median per-tile kernel sum (flux scale)")
        conv_hdr["COLCLEAN"] = (col_clean, "column medians removed, noisy columns zeroed")
        conv_hdr["NBADCOL"] = (len(badcols), "columns zeroed by clean_columns")
        conv_hdr["BADCOLS"] = (tc.badcols_card(badcols), "zeroed cols")
        add_member(tar, f"conv_{f['name']}", tc.fits_bytes(diff.astype(np.float32), conv_hdr))

        interp_hdr = hdr.copy()
        interp_hdr["INTERPMD"] = (interp_mode, "full image, or header-only for WCS")
        interp_data = f["interp"] if interp_mode == "full" else np.zeros((1, 1), np.float32)
        add_member(tar, f"interp_{f['name']}", tc.fits_bytes(interp_data, interp_hdr))

        sky_hdr = fits.Header()
        for key in ("STARTTJD", "MIDTJD", "ENDTJD", "QUAL_BIT"):
            if key in f["hdr"]:
                sky_hdr[key] = f["hdr"][key]
        sky_hdr["SKYBOX"] = (sky_box, "block size of this sky-level mesh [pix]")
        sky_hdr["SCIX0"] = (tc.SCI_X.start, "0-based x of first science pixel")
        sky_hdr["SCIY0"] = (tc.SCI_Y.start, "0-based y of first science pixel")
        diff_hdr = fits.Header()
        diff_hdr["DIFFBOX"] = (diff_box, "Background2D box of this mesh [pix]")
        buf = io.BytesIO()
        fits.HDUList([fits.PrimaryHDU(f["sky_mesh"], sky_hdr),
                      fits.ImageHDU(diff_mesh, diff_hdr, name="DIFFBKG")]).writeto(buf)
        add_member(tar, f"bkg_{f['name']}", buf.getvalue())

        stats.append((f["name"], f["dx"], f["dy"], f["bkg_med"], diff_rms, scale, resid,
                      f["raw_bytes"], f["stored_bytes"]))
os.replace(tar_tmp, os.path.join(slice_dir, "images.tar"))

with open(os.path.join(slice_dir, "dates"), "w") as out:
    for f in frames:
        out.write(f"{f['name']} {f['tjd']:.8f}\n")
with open(os.path.join(slice_dir, ".quality"), "w") as out:
    for f in frames:
        out.write(f"{f['name']} {f['qual']} {f['tjd']:.8f}\n")
with open(os.path.join(slice_dir, ".stats.tsv"), "w") as out:
    out.write("name\tshift_x\tshift_y\tbkg_med\tdiff_rms\tkscale\tconv_std\traw_bytes\tstored_bytes\n")
    for row in stats:
        out.write(f"{row[0]}\t{row[1]:.4f}\t{row[2]:.4f}\t{row[3]:.3f}\t{row[4]:.3f}\t"
                  f"{row[5]:.5f}\t{row[6]:.4f}\t{row[7]}\t{row[8]}\n")

if keep_raw:
    if q is not None:  # keep only the compressed copies
        for f in frames:
            path = os.path.join(cache, f["name"])
            if os.path.exists(path):
                os.remove(path)
else:
    shutil.rmtree(cache, ignore_errors=True)

open(done_path, "w").close()
print(f"Done {args.orbit}/{args.slice}: {len(frames)} frames in {time.time() - t0:.1f}s")
