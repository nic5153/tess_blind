"""Shared helpers for the local TICA upstream processing.

These scripts rebuild, on a local machine, the products that the HPCC
image-subtraction step used to leave under $DATA_DIR:

    $DATA_DIR/s00XX/camC-ccdD/ref.fits
    $DATA_DIR/s00XX/camC-ccdD/oXX/rms_list_tica
    $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/dates
    $DATA_DIR/s00XX/camC-ccdD/oXX/sliceNNNN/images.tar   (conv_, interp_, bkg_<name>)

so the existing make_rms / source / analysis scripts only need path changes.
"""
import io
import os
import re
import time
import urllib.error
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip
from astropy.wcs import WCS
from photutils.background import Background2D, MedianBackground

TICA_URL = os.environ.get("TICA_URL", "https://archive.stsci.edu/hlsps/tica")
NAME_RE = re.compile(
    r"hlsp_tica_tess_ffi_s(\d{4})-(o[12])-(\d{8})-cam(\d)-ccd(\d)_tess_v01_img\.fits"
)

# Science pixels from the TICA SCIPIXS keyword '[45:2092,1:2048]' (1-based, inclusive).
SCI_Y = slice(0, 2048)
SCI_X = slice(44, 2092)

WCS_KEY_RE = re.compile(
    r"^(WCSAXES|CTYPE\d|CRPIX\d|CRVAL\d|CDELT\d|CUNIT\d|CD\d_\d|PC\d_\d|"
    r"LONPOLE|LATPOLE|RADESYS|(A|B|AP|BP)_(\d+_\d+|ORDER))$"
)


def env_str(name, default):
    value = os.environ.get(name, "")
    return value if value != "" else default


def env_int(name, default):
    try:
        return int(env_str(name, default))
    except ValueError:
        return default


def env_float(name, default):
    try:
        return float(env_str(name, default))
    except ValueError:
        return default


def compress_q():
    """fpack-style quantization level for the lossy-compression test, or None."""
    value = os.environ.get("COMPRESS_Q", "").strip()
    return float(value) if value else None


def data_dir():
    return env_str("DATA_DIR", "/mnt/d/tess_blind_data/tica")


def raw_cache():
    return env_str("RAW_CACHE", "/mnt/d/tess_blind_data/raw_cache")


def ccd_dir(sector, cam, ccd):
    return os.path.join(data_dir(), f"s{sector:04d}", f"cam{cam}-ccd{ccd}")


def frame_url(sector, cam, ccd, name):
    return f"{TICA_URL}/s{sector:04d}/cam{cam}-ccd{ccd}/{name}"


def manifest_path(sector, cam, ccd):
    return os.path.join(ccd_dir(sector, cam, ccd), "manifest.tsv")


def read_manifest(sector, cam, ccd):
    path = manifest_path(sector, cam, ccd)
    if not os.path.exists(path):
        raise FileNotFoundError(f"missing manifest {path}; run tica_manifest.py first")
    rows = []
    with open(path) as f:
        header = f.readline().split()
        for line in f:
            parts = line.split()
            if len(parts) != len(header):
                continue
            row = dict(zip(header, parts))
            row["cadence"] = int(row["cadence"])
            rows.append(row)
    return rows


def download(url, dest, retries=4, timeout=180):
    """Download url to dest atomically (via dest.part), retrying on failure."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    tmp = dest + ".part"
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp, open(tmp, "wb") as out:
                while True:
                    chunk = resp.read(1 << 20)
                    if not chunk:
                        break
                    out.write(chunk)
            os.replace(tmp, dest)
            return dest
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            if os.path.exists(tmp):
                os.remove(tmp)
            if attempt == retries:
                raise RuntimeError(f"download failed after {retries} tries: {url} ({exc})")
            time.sleep(5 * attempt)


def pool_dir(sector, cam, ccd):
    """Where tica_prefetch.py stores frames ahead of time (for compute nodes
    without internet access)."""
    return os.path.join(raw_cache(), f"s{sector:04d}", f"cam{cam}-ccd{ccd}", "pool")


def fetch_frame(sector, cam, ccd, name, cache_dir, take=True):
    """Local path of a frame in cache_dir: taken from the prefetch pool if there
    (moved when take=True, copied otherwise), else downloaded from MAST."""
    dest = os.path.join(cache_dir, name)
    if os.path.exists(dest):
        return dest
    pooled = os.path.join(pool_dir(sector, cam, ccd), name)
    if os.path.exists(pooled):
        import shutil

        os.makedirs(cache_dir, exist_ok=True)
        if take:
            shutil.move(pooled, dest)
        else:
            shutil.copyfile(pooled, dest)
        return dest
    download(frame_url(sector, cam, ccd, name), dest)
    return dest


def rice_roundtrip(data, header, q, keep_path=None):
    """Lossy RICE compression exactly as fpack does it (quantize level q,
    SUBTRACTIVE_DITHER_1), then decompress. Returns (data, compressed_bytes)."""
    hdu = fits.CompImageHDU(
        data=data,
        header=header,
        compression_type="RICE_1",
        quantize_level=q,
        quantize_method=1,
        dither_seed=-1,  # seed from the data checksum: reproducible
    )
    buf = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(buf)
    nbytes = buf.tell()
    if keep_path is not None:
        with open(keep_path, "wb") as f:
            f.write(buf.getvalue())
    buf.seek(0)
    with fits.open(buf) as hdul:
        out = np.array(hdul[1].data, dtype=np.float32)
    return out, nbytes


def load_frame(path, q=None, keep_compressed=None):
    """Read a TICA frame as float32. With q set, pass it through fpack-style
    lossy compression first. Returns (data, header, raw_bytes, stored_bytes)."""
    with fits.open(path, memmap=False) as hdul:
        data = np.array(hdul[0].data, dtype=np.float32)
        header = hdul[0].header.copy()
    raw_bytes = os.path.getsize(path)
    stored_bytes = raw_bytes
    if q is not None:
        data, stored_bytes = rice_roundtrip(data, header, q, keep_compressed)
    return data, header, raw_bytes, stored_bytes


def frame_wcs(header):
    import warnings

    from astropy.wcs import FITSFixedWarning

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        return WCS(header, relax=True)


def lstsq_kernel(image, ref, xs, ys, half=2, stamp=6, ksize=13, nclip=3):
    """Real-space PSF-matching kernel (Bramich 2008 delta-function basis), the
    approach ISIS-style difference imaging uses: solve image ~ ref (x) k + c by
    least squares over stamps around stars, with iterative 4-sigma clipping of
    variable/saturated pixels. The kernel sum carries the flux scale.
    Returns a ksize x ksize kernel, or None if there are too few stars."""
    n = 2 * half + 1
    rows, targets = [], []
    for x, y in zip(xs, ys):
        ix, iy = int(round(x)), int(round(y))
        tgt = image[iy - stamp:iy + stamp + 1, ix - stamp:ix + stamp + 1]
        if tgt.shape != (2 * stamp + 1, 2 * stamp + 1):
            continue
        cols = [ref[iy - stamp - v:iy + stamp + 1 - v, ix - stamp - u:ix + stamp + 1 - u].ravel()
                for v in range(-half, half + 1) for u in range(-half, half + 1)]
        cols.append(np.ones(tgt.size))
        rows.append(np.array(cols).T)
        targets.append(tgt.ravel())
    if len(rows) < 5:
        return None
    A = np.vstack(rows).astype(np.float64)
    b = np.concatenate(targets).astype(np.float64)
    good = np.isfinite(b) & np.all(np.isfinite(A), axis=1)
    coef = None
    for _ in range(nclip):
        coef, *_ = np.linalg.lstsq(A[good], b[good], rcond=None)
        resid = b - A @ coef
        sig = 1.4826 * np.median(np.abs(resid[good]))
        good = np.isfinite(b) & (np.abs(resid) < 4 * sig)
    kern = np.zeros((ksize, ksize))
    c = ksize // 2
    # column (u, v) holds ref(y - v, x - u): a true convolution with kern[v, u]
    kern[c - half:c + half + 1, c - half:c + half + 1] = coef[:n * n].reshape(n, n)
    return kern


def wcs_shift(ref_wcs, wcs, ngrid=7):
    """Median (dx, dy) such that a source at ref pixel p sits at frame pixel p + d."""
    xs = np.linspace(SCI_X.start + 100, SCI_X.stop - 100, ngrid)
    ys = np.linspace(SCI_Y.start + 100, SCI_Y.stop - 100, ngrid)
    gx, gy = np.meshgrid(xs, ys)
    ra, dec = ref_wcs.all_pix2world(gx.ravel(), gy.ravel(), 0)
    fx, fy = wcs.all_world2pix(ra, dec, 0)
    return float(np.median(fx - gx.ravel())), float(np.median(fy - gy.ravel()))


def register(data, dx, dy):
    """Resample a frame onto the reference grid given its wcs_shift (dx, dy).
    Shifts below MIN_SHIFT_PIX are skipped (TESS jitter is usually tiny)."""
    from scipy import ndimage

    if max(abs(dx), abs(dy)) < env_float("MIN_SHIFT_PIX", 0.02):
        return data
    return ndimage.shift(data, (-dy, -dx), order=3, mode="nearest").astype(np.float32)


def saturation_level(header):
    """Conservative per-pixel saturation in a TICA co-add: ~1.5e5 e- per 2 s
    exposure (TESS full well ~2e5 e-) times the number of 2 s exposures."""
    exptime = float(header.get("EXPTIME", 475.2))
    return env_float("SAT_E_PER_2S", 1.5e5) * exptime / 2.0


def science_background(data, box=None, return_mesh=False):
    """Background2D over the science region; returns a full-frame background array
    (zero outside the science pixels), its median and median rms, and with
    return_mesh=True the low-resolution background mesh as well."""
    box = box or env_int("BKG_BOX", 64)
    bkg = Background2D(
        data[SCI_Y, SCI_X],
        (box, box),
        filter_size=(3, 3),
        sigma_clip=SigmaClip(sigma=3.0),
        bkg_estimator=MedianBackground(),
    )
    full = np.zeros_like(data, dtype=np.float32)
    full[SCI_Y, SCI_X] = bkg.background
    out = (full, float(np.median(bkg.background)), float(np.median(bkg.background_rms)))
    if return_mesh:
        out += (np.asarray(bkg.background_mesh, dtype=np.float32),)
    return out


def clean_columns(diff, nsig=None, grow=None):
    """Remove column artifacts from a difference image, in place:
      1) subtract each science column's median (smear-like column offsets);
      2) zero columns whose robust noise exceeds nsig x the typical column's
         (bleed/smear columns of bright stars: full-height lines in the daily
         RMS images that otherwise become false candidates), plus `grow`
         neighbors on each side.
    Returns the list of zeroed full-frame column indices."""
    nsig = env_float("COLMASK_NSIG", 5.0) if nsig is None else nsig
    grow = env_int("COLMASK_GROW", 1) if grow is None else grow
    sci = diff[SCI_Y, SCI_X]
    sci -= np.median(sci, axis=0)[None, :]
    noise = 1.4826 * np.median(np.abs(sci), axis=0)
    bad = noise > nsig * np.median(noise)
    if grow > 0 and bad.any():
        idx = np.where(bad)[0]
        for g in range(1, grow + 1):
            bad[np.clip(idx - g, 0, bad.size - 1)] = True
            bad[np.clip(idx + g, 0, bad.size - 1)] = True
    sci[:, bad] = 0.0
    diff[SCI_Y, SCI_X] = sci
    return (np.where(bad)[0] + SCI_X.start).tolist()


def badcols_card(cols, maxlen=40):
    """Compact 'a-b,c' string of column indices for a FITS header card."""
    if not cols:
        return ""
    runs, start, prev = [], cols[0], cols[0]
    for c in cols[1:] + [None]:
        if c is not None and c == prev + 1:
            prev = c
            continue
        runs.append(f"{start}" if start == prev else f"{start}-{prev}")
        if c is not None:
            start = prev = c
    text = ",".join(runs)
    return text if len(text) <= maxlen else text[:maxlen - 3] + "..."


def block_median(data, box):
    """Median of each box x box block of the science region (a cheap sky map)."""
    sci = data[SCI_Y, SCI_X]
    ny, nx = sci.shape[0] // box, sci.shape[1] // box
    blocks = sci[:ny * box, :nx * box].reshape(ny, box, nx, box).transpose(0, 2, 1, 3)
    return np.median(blocks.reshape(ny, nx, box * box), axis=2).astype(np.float32)


def mesh_to_image(mesh, box, shape):
    """Expand a science-region mesh (box x box blocks) back to a full-frame image."""
    from scipy import ndimage

    full = np.zeros(shape, dtype=np.float32)
    up = ndimage.zoom(mesh, box, order=1, mode="nearest", grid_mode=True)
    ny = min(up.shape[0], SCI_Y.stop - SCI_Y.start)
    nx = min(up.shape[1], SCI_X.stop - SCI_X.start)
    full[SCI_Y.start:SCI_Y.start + ny, SCI_X.start:SCI_X.start + nx] = up[:ny, :nx]
    return full


def tile_index(x, y, ntile):
    """Kernel/PSF tile (ty * ntile + tx) containing 0-based pixel (x, y)."""
    xedges = np.linspace(SCI_X.start, SCI_X.stop, ntile + 1)
    yedges = np.linspace(SCI_Y.start, SCI_Y.stop, ntile + 1)
    tx = np.clip(np.searchsorted(xedges, x, side="right") - 1, 0, ntile - 1)
    ty = np.clip(np.searchsorted(yedges, y, side="right") - 1, 0, ntile - 1)
    return ty * ntile + tx


def replace_wcs(header, wcs_header):
    """Swap the WCS in header for wcs_header (used once a frame is resampled
    onto the reference pixel grid)."""
    out = header.copy()
    for key in list(out.keys()):
        if WCS_KEY_RE.match(key):
            del out[key]
    out.update(wcs_header)
    return out


def tiles(ntile):
    """Yield (ty, tx, yslice, xslice) tiles covering the science region."""
    ny = SCI_Y.stop - SCI_Y.start
    nx = SCI_X.stop - SCI_X.start
    yedges = np.linspace(SCI_Y.start, SCI_Y.stop, ntile + 1).astype(int)
    xedges = np.linspace(SCI_X.start, SCI_X.stop, ntile + 1).astype(int)
    assert ny > 0 and nx > 0
    for ty in range(ntile):
        for tx in range(ntile):
            yield ty, tx, slice(yedges[ty], yedges[ty + 1]), slice(xedges[tx], xedges[tx + 1])


def stamp_stack(image, xs, ys, size):
    """Flux-normalized median of size x size stamps centered on the nearest pixel
    to each (x, y). Returns (psf or None, per-star stamp sums)."""
    half = size // 2
    stamps, sums = [], []
    for x, y in zip(xs, ys):
        ix, iy = int(round(x)), int(round(y))
        stamp = image[iy - half : iy + half + 1, ix - half : ix + half + 1]
        if stamp.shape != (size, size) or not np.all(np.isfinite(stamp)):
            sums.append(np.nan)
            continue
        total = float(stamp.sum())
        sums.append(total)
        if total > 0:
            stamps.append(stamp / total)
    if len(stamps) < 3:
        return None, np.array(sums)
    psf = np.median(np.array(stamps), axis=0)
    psf = np.clip(psf, 0, None)
    psf /= psf.sum()
    return psf, np.array(sums)


def fits_bytes(data, header):
    buf = io.BytesIO()
    fits.PrimaryHDU(data=data, header=header).writeto(buf)
    return buf.getvalue()
