import argparse
import io
import os
import re
import tarfile
import traceback
import warnings
from io import StringIO
from multiprocessing import Pool

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as astrou
from astropy.io import fits
from astropy.time import Time
from astropy.wcs import WCS
from matplotlib.gridspec import GridSpec
from scipy.optimize import curve_fit
from tqdm import tqdm
from tess_stars2px import tess_stars2px_reverse_function_entry


font = {"weight": "bold", "size": 6}
mpl.rc("font", **font)

parser = argparse.ArgumentParser()
parser.add_argument("--photfile", required=True, type=str, help="phot.data file")
parser.add_argument("--singleproc", action="store_true")
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--N", type=int, default=4)
parser.add_argument("--rmsthresh", type=float, default=100.0)
parser.add_argument("--data-dir", default=os.environ.get("DATA_DIR", "/lustre/research/mfausnau/data/tica"))
args = parser.parse_args()

cam = int(args.photfile.replace("/phot.data", "").split("/")[-1][3])
ccd = int(args.photfile.replace("/phot.data", "").split("/")[-1][-1])
print(cam, ccd)

sourcedata = np.loadtxt(
    args.photfile,
    dtype=[
        ("fcol", "f8"),
        ("frow", "f8"),
        ("icol", "i8"),
        ("irow", "i8"),
        ("fname", "<U256"),
        ("flag", "i4"),
    ],
)
if sourcedata.size == 1:
    sourcedata = sourcedata[np.newaxis]


def load_dates():
    alldates = ""
    pattern = f"{args.data_dir}/s{args.sector:04}/cam{cam}-ccd{ccd}/o??/slice*/dates"
    for datesfile in glob_sorted(pattern):
        thisorbit = datesfile.split("/")[-3]
        with open(datesfile, "r") as f:
            for line in f:
                alldates += f"{thisorbit} {line}"

    if not alldates.strip():
        return np.array([], dtype=[("detorbit", "<U8"), ("name", "<U128"), ("dates", "f8")])

    dates = np.unique(
        np.genfromtxt(
            StringIO(alldates),
            dtype=[("detorbit", "<U8"), ("name", "<U128"), ("dates", "f8")],
        )
    )
    if dates.size == 1:
        dates = dates[np.newaxis]
    dates.sort(order="dates")
    return dates


def glob_sorted(pattern):
    import glob

    return sorted(glob.glob(pattern))


dates = load_dates()

# _cleaned files written by clean_handmade_lc.py:
CLEAN_BTJD = 0
CLEAN_TJD = 1
CLEAN_CTS_PER_S = 2
CLEAN_E_CTS_PER_S = 3
CLEAN_MAG = 4
CLEAN_E_MAG = 5
CLEAN_BKG = 6
CLEAN_BKG_MODEL = 7
CLEAN_BKG2 = 8
CLEAN_E_BKG2 = 9

# Raw phot2 light curves used by this no-cleaned plotter.
RAW_TJD = 0
RAW_DIFF_COUNTS = 1
RAW_E_DIFF_COUNTS = 2
RAW_BKG_COUNTS = 6


def ctstomag(cts):
    with np.errstate(divide="ignore", invalid="ignore"):
        return -2.5 * np.log10(cts) + 20.44


def gaussian_nobackground(xy, amplitude, xo, yo, sigma_x, sigma_y, theta):
    x, y = xy
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta) ** 2) / (2 * sigma_x**2) + (np.sin(theta) ** 2) / (2 * sigma_y**2)
    b = -(np.sin(2 * theta)) / (4 * sigma_x**2) + (np.sin(2 * theta)) / (4 * sigma_y**2)
    c = (np.sin(theta) ** 2) / (2 * sigma_x**2) + (np.cos(theta) ** 2) / (2 * sigma_y**2)
    g = amplitude * np.exp(-(a * ((x - xo) ** 2) + 2 * b * (x - xo) * (y - yo) + c * ((y - yo) ** 2)))
    return g.ravel()


def parse_candidate_day(lcfile):
    match = re.search(r"outcatrms_day_(\d+)", os.path.basename(lcfile))
    if match is None:
        return None
    return int(match.group(1))


def raw_light_curve(lcfile):
    d = np.genfromtxt(lcfile)
    if d.ndim == 1:
        d = d[np.newaxis, :]
    d = d[(~np.isnan(d[:, RAW_TJD])) & (~np.isnan(d[:, RAW_DIFF_COUNTS])) & (~np.isnan(d[:, RAW_E_DIFF_COUNTS]))]
    d = d[np.argsort(d[:, RAW_TJD])]

    if len(d) < 2:
        raise ValueError(f"Not enough points in {lcfile}")

    diffs = np.diff(d[:, RAW_TJD])
    diffs = diffs[np.isfinite(diffs) & (diffs > 0)]
    if len(diffs) == 0:
        tdiff = 200.0
    else:
        tdiff = float(np.nanmedian(np.round(diffs * 24 * 60 * 60)))

    scale = tdiff * 0.8 * 0.99

    # Legacy cleaned light curves used columns 0,2,3,4,5,6 directly as
    # BTJD/cts_per_s/e_cts_per_s/mag/e_mag/bkg. The no-cleaned files are raw
    # phot2 outputs, so we keep the same cts/s scale by dividing counts by the
    # exposure factor used in clean_handmade_lc.py.
    xvals = d[:, RAW_TJD]
    yvals = -d[:, RAW_DIFF_COUNTS] / scale
    yerr = d[:, RAW_E_DIFF_COUNTS] / scale
    ybkg = d[:, RAW_BKG_COUNTS] / scale if d.shape[1] > RAW_BKG_COUNTS else np.full_like(xvals, np.nan)
    return xvals, yvals, yerr, ybkg, tdiff, scale


def background_lc_path(lcfile):
    parts = lcfile.replace("\\", "/").split("/")
    try:
        idx = parts.index("lc")
    except ValueError:
        return lcfile.replace("/lc/", "/bkg_phot/lc/")
    return "/".join(parts[:idx] + ["bkg_phot", "lc"] + parts[idx + 1 :])


def align_by_time(source_times, value_times, values):
    aligned = np.full_like(source_times, np.nan, dtype=float)
    source_key = np.round(source_times, 8)
    value_key = np.round(value_times, 8)
    _, source_idx, value_idx = np.intersect1d(source_key, value_key, return_indices=True)
    aligned[source_idx] = values[value_idx]
    return aligned


def find_source_row(lcfile):
    srcname = os.path.basename(lcfile)
    candidates = [
        f"lc/{srcname}",
        srcname,
    ]
    names = sourcedata["fname"].astype(str)
    for candidate in candidates:
        match = names == candidate
        if np.any(match):
            return sourcedata[match][0]
    raise ValueError(f"No phot.data row for {srcname}")


def choose_peak_time(xvals, yvals, candidate_day):
    if candidate_day is not None:
        daymask = (xvals >= candidate_day) & (xvals < candidate_day + 1)
        if np.any(daymask):
            subx = xvals[daymask]
            suby = yvals[daymask]
            return subx[np.nanargmax(suby)]
    return xvals[np.nanargmax(yvals)]


def nearest_date(tref):
    if dates.size == 0:
        return None
    idx = np.nanargmin(np.abs(dates["dates"] - tref))
    return dates[idx]


def name_in_dates(dates_path, target_name):
    with open(dates_path, "r") as f:
        for line in f:
            parts = line.split()
            if parts and parts[0] == target_name:
                return True
    return False


def fits_from_tar(tar_path, target_names):
    out = {}
    with tarfile.open(tar_path, "r") as tar:
        for member in tar.getmembers():
            basename = os.path.basename(member.name)
            if basename not in target_names:
                continue
            fobj = tar.extractfile(member)
            if fobj is None:
                continue
            out[basename] = fobj.read()
            if len(out) == len(target_names):
                break
    return out


def load_detection_images(tref, scale):
    det = nearest_date(tref)
    if det is None:
        raise RuntimeError("No dates files were loaded")

    base_path = f"{args.data_dir}/s{args.sector:04}/cam{cam}-ccd{ccd}/{det['detorbit']}"
    target_name = str(det["name"])
    conv_name = f"conv_{target_name}"
    interp_name = f"interp_{target_name}"

    for dates_file in glob_sorted(f"{base_path}/slice*/dates"):
        if not name_in_dates(dates_file, target_name):
            continue
        slice_dir = os.path.dirname(dates_file)

        loose_conv = os.path.join(slice_dir, conv_name)
        loose_interp = os.path.join(slice_dir, interp_name)
        if os.path.exists(loose_conv) and os.path.exists(loose_interp):
            with fits.open(loose_interp) as hdul:
                header = hdul[0].header
            with fits.open(loose_conv) as hdul:
                data = -hdul[0].data / scale
            return data, header, det

        tar_path = os.path.join(slice_dir, "images.tar")
        if not os.path.exists(tar_path):
            continue
        packed = fits_from_tar(tar_path, {conv_name, interp_name})
        if conv_name not in packed or interp_name not in packed:
            continue
        with fits.open(io.BytesIO(packed[interp_name])) as hdul:
            header = hdul[0].header
        with fits.open(io.BytesIO(packed[conv_name])) as hdul:
            data = -hdul[0].data / scale
        return data, header, det

    raise RuntimeError(f"Could not find conv/interp image for {target_name}")


def local_camccd_dir():
    return os.path.dirname(args.photfile.rstrip("/"))


def daily_rms_files(candidate_day):
    root = local_camccd_dir()
    if candidate_day is None:
        pattern = os.path.join(root, "o??", "rms_day_*.fits")
    else:
        pattern = os.path.join(root, "o??", f"rms_day_{candidate_day}.fits")
    return glob_sorted(pattern)


def all_daily_rms_files():
    return glob_sorted(os.path.join(local_camccd_dir(), "o??", "rms_day_*.fits"))


def display_window(ax, data, icol, irow, fcol, frow, imsize, vmin=None, vmax=None, title=None):
    sub = data[(irow - imsize) : (irow + imsize + 1), (icol - imsize) : (icol + imsize + 1)]
    if vmin is None or vmax is None:
        flat = np.ravel(sub[np.isfinite(sub)])
        if len(flat):
            vmin = np.nanpercentile(flat, 5)
            vmax = np.nanpercentile(flat, 95)
        else:
            vmin, vmax = 0, 1
    ax.imshow(sub, cmap="Greys", origin="lower", vmin=vmin, vmax=vmax)
    circle = plt.Circle(((fcol - icol) + imsize, (frow - irow) + imsize), 2.5, fill=False)
    ax.add_artist(circle)
    if title:
        ax.set_title(title)
    ax.axis("off")


def coords_from_image_position(header, col, row):
    try:
        src_sc = WCS(header).pixel_to_world(col, row)
        ra = src_sc.ra.deg
        dec = src_sc.dec.deg
    except Exception:
        ra, dec, _ = tess_stars2px_reverse_function_entry(args.sector, cam, ccd, col, row)
        src_sc = SkyCoord(ra, dec, unit="deg")
    return ra, dec, src_sc.transform_to("galactic")


def fit_gaussian_panel(fig, gs, image, header, src, imsize):
    frow = src["frow"]
    fcol = src["fcol"]
    icol = int(round(src["icol"]))
    irow = int(round(src["irow"]))

    axs = [
        fig.add_subplot(gs[0, -2]),
        fig.add_subplot(gs[0, -1]),
        fig.add_subplot(gs[1, -2]),
        fig.add_subplot(gs[1, -1]),
    ]
    axs[0].set_title("data")
    subdata = np.copy(image[(irow - imsize) : (irow + imsize + 1), (icol - imsize) : (icol + imsize + 1)])

    constraint = 4
    imdata = np.copy(image[(irow - constraint) : (irow + constraint + 1), (icol - constraint) : (icol + constraint + 1)])
    if imdata.shape[0] != imdata.shape[1]:
        raise RuntimeError("candidate is too close to the edge")

    x0 = (frow - irow) + constraint
    y0 = (fcol - icol) + constraint
    f0 = imdata[int(round(x0)), int(round(y0))]
    bsizepix = 1.2
    # Legacy bounds used f0 directly. Raw difference images can put a negative
    # or non-finite value at the starting pixel, which inverts the amplitude
    # bounds and makes scipy fail before the original ellipticity test can run.
    if not np.isfinite(f0) or f0 <= 0:
        positive_pixels = imdata[np.isfinite(imdata) & (imdata > 0)]
        f0 = np.nanmax(positive_pixels) if positive_pixels.size else 1.0
    minf0 = max(0, f0 * 0.5)
    maxf0 = min(1e6, f0 * 1.5)
    if maxf0 <= minf0:
        minf0 = 0
        maxf0 = 1e9
    bounds = [
        (minf0, x0 - constraint, y0 - constraint, bsizepix * 0.5, bsizepix * 0.5, 0),
        (maxf0, x0 + constraint, y0 + constraint, bsizepix * 5, bsizepix * 5, 2 * np.pi),
    ]

    fitdata = np.transpose(imdata)
    xarr = np.arange(2 * constraint + 1)
    yarr = np.arange(2 * constraint + 1)
    subx, suby = np.meshgrid(xarr, yarr)
    try:
        popt, _ = curve_fit(
            gaussian_nobackground,
            (subx, suby),
            fitdata.ravel(),
            p0=[f0, x0, y0, bsizepix, bsizepix, 0],
            bounds=bounds,
        )
    except RuntimeError:
        vmax = np.nanmax(imdata) if np.isfinite(imdata).any() else np.nanmax(subdata)
        vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1
        display_window(axs[0], image, icol, irow, fcol, frow, imsize, vmin=0, vmax=vmax, title="data")
        axs[1].set_title("model")
        axs[1].text(0.05, 0.5, "fit failed", transform=axs[1].transAxes)
        axs[1].axis("off")
        axs[2].set_title("residuals")
        axs[2].text(0.05, 0.5, "fit failed", transform=axs[2].transAxes)
        axs[2].axis("off")
        axs[3].text(0, 0.7, "Gaussian fit values")
        axs[3].text(0, 0.5, "fit failed")
        axs[3].text(0, 0.1, "no convergence")
        axs[3].axis("off")
        ra, dec, src_gal = coords_from_image_position(header, fcol, frow)
        return ra, dec, src_gal, np.nan

    true_row = popt[1] + (irow - constraint)
    true_col = popt[2] + (icol - constraint)
    ra, dec, src_gal = coords_from_image_position(header, true_col, true_row)

    model = gaussian_nobackground((subx, suby), *popt).reshape(fitdata.shape)
    vmax = np.nanmax(model) if np.isfinite(model).any() else np.nanmax(subdata)
    vmax = vmax if np.isfinite(vmax) and vmax > 0 else 1

    display_window(axs[0], image, icol, irow, fcol, frow, imsize, vmin=0, vmax=vmax, title="data")
    axs[1].imshow(model, cmap="Greys", origin="lower", vmin=0, vmax=vmax)
    axs[1].set_title("model")
    axs[1].axis("off")
    axs[2].imshow(np.abs(fitdata - model), cmap="Greys", origin="lower", vmin=0, vmax=vmax)
    axs[2].set_title("residuals")
    axs[2].axis("off")
    axs[3].text(0, 0.7, "Gaussian fit values")
    axs[3].text(0, 0.5, "Amp, col, row, rad1, rad2, theta")
    axs[3].text(0, 0.1, f"{popt[0]:5.4e},{true_col:6.1f},{true_row:6.1f},\n{popt[3]:5.4e},{popt[4]:5.4e},{popt[5]:5.4e}")
    axs[3].axis("off")

    ellip_ratio = popt[-3] / popt[-2]
    return ra, dec, src_gal, ellip_ratio


def add_reference_panel(fig, gs, scale, src, imsize):
    rootfolder = f"{args.data_dir}/s{args.sector:04}/cam{cam}-ccd{ccd}"
    reffile = f"{rootfolder}/ref.fits"
    ax = fig.add_subplot(gs[0, 2])
    try:
        with fits.open(reffile) as hdul:
            data = hdul[0].data / scale
        display_window(ax, data, int(round(src["icol"])), int(round(src["irow"])), src["fcol"], src["frow"], imsize, vmin=0, vmax=400, title="ref")
    except Exception as e:
        ax.text(0, 0.5, f"ref unavailable\n{e}")
        ax.axis("off")


def add_daily_rms_panels(fig, gs, candidate_day, scale, src, imsize, detorbit_hint=None):
    frow = src["frow"]
    fcol = src["fcol"]
    icol = int(round(src["icol"]))
    irow = int(round(src["irow"]))

    allrms = all_daily_rms_files()
    detfiles = daily_rms_files(candidate_day)
    if not allrms or not detfiles:
        warnings.warn(f"Could not find daily RMS image for day {candidate_day}")
        return None

    # Legacy daily RMS selection used the first rms_day match. If the image
    # dates lookup succeeded, prefer the same orbit while preserving the old
    # first-match fallback.
    if detorbit_hint is not None:
        same_detection_orbit = [f for f in detfiles if os.path.basename(os.path.dirname(f)) == detorbit_hint]
        detfile = same_detection_orbit[0] if same_detection_orbit else detfiles[0]
    else:
        detfile = detfiles[0]
    detorbit = os.path.basename(os.path.dirname(detfile))
    same_orbit = [f for f in allrms if os.path.basename(os.path.dirname(f)) == detorbit]
    context = same_orbit if len(same_orbit) > 1 else allrms
    context = sorted(context)

    detidx = context.index(detfile) if detfile in context else 0
    start = max(0, detidx - 1)
    panel_files = context[start : start + 4]
    if detfile not in panel_files:
        panel_files = [detfile] + [f for f in panel_files if f != detfile]
    panel_files = panel_files[:4]

    axs = []
    for i, _ in enumerate(panel_files):
        rowno = i % 2
        colno = int((i / 2) % 2)
        axs.append(fig.add_subplot(gs[colno, rowno]))

    images = []
    detim = None
    otherim = []
    for ffile in panel_files:
        with fits.open(ffile) as hdul:
            data = hdul[0].data / scale
        images.append((ffile, data))
        if ffile == detfile:
            detim = np.copy(data)
        else:
            otherim.append(data)

    if detim is None:
        return None

    for ax, (ffile, data) in zip(axs, images):
        display_window(ax, data, icol, irow, fcol, frow, imsize, vmin=0, vmax=5, title=os.path.basename(ffile).replace(".fits", ""))

    if otherim:
        median_other = np.nanmedian(np.array(otherim), axis=0)
        detsub = detim - median_other
        ax = fig.add_subplot(gs[1, 2])
        display_window(ax, detsub, icol, irow, fcol, frow, imsize, vmin=0, vmax=5, title="daily RMS - median")
        return detsub

    return detim


def plot(lcfile):
    try:
        if "_cleaned" in lcfile or lcfile.endswith(".png"):
            return

        candidate_day = parse_candidate_day(lcfile)
        src = find_source_row(lcfile)
        frow = src["frow"]
        fcol = src["fcol"]
        icol = int(round(src["icol"]))
        irow = int(round(src["irow"]))

        xvals, yvals, yerr, ybkg, tdiff, scale = raw_light_curve(lcfile)
        bkgfile = background_lc_path(lcfile)
        if os.path.exists(bkgfile):
            # Legacy assumption after matchbkgtolc.py:
            # bxvals, byvals, byerr, bybkg, _, _ = raw_light_curve(bkgfile)
            #
            # Keep the same background product, but align by rounded TJD so
            # analysis.sbatch can run even when the external matcher is absent
            # or leaves extra background cadences.
            bxvals_raw, byvals_raw, byerr_raw, bybkg_raw, _, _ = raw_light_curve(bkgfile)
            bxvals = xvals
            byvals = align_by_time(xvals, bxvals_raw, byvals_raw)
            byerr = align_by_time(xvals, bxvals_raw, byerr_raw)
            bybkg = align_by_time(xvals, bxvals_raw, bybkg_raw)
        else:
            bxvals, byvals, byerr, bybkg = xvals, np.full_like(xvals, np.nan), np.full_like(xvals, np.nan), np.full_like(xvals, np.nan)

        tpeak = choose_peak_time(xvals, yvals, candidate_day)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ra, dec, _ = tess_stars2px_reverse_function_entry(args.sector, cam, ccd, fcol, frow)
        src_gal = SkyCoord(ra, dec, unit="deg").transform_to("galactic")

        flaggedlist = []
        rmslen = 20
        runrms = []
        rms_signal = bybkg if np.isfinite(bybkg).any() else ybkg
        # Legacy flagging threshold used the background signal directly:
        # if currms > 100: flaggedlist.extend(curxvals)
        #
        # In raw no-cleaned files column 6 includes a large absolute background
        # pedestal. Subtracting the local median preserves the same RMS test but
        # applies it to background variation instead of the DC level.
        if np.any(np.isfinite(rms_signal)):
            rms_signal = rms_signal - np.nanmedian(rms_signal)
        for i in range(max(0, len(xvals) - rmslen)):
            curvals = rms_signal[i : i + rmslen]
            curxvals = xvals[i : i + rmslen]
            if len(curxvals) == 0 or (curxvals.max() - curxvals.min()) > 1:
                currms = np.nan
            else:
                currms = np.sqrt(np.nanmean(curvals**2))
            runrms.append((np.nanmean(curxvals), currms))
            if np.isfinite(currms) and currms > args.rmsthresh:
                flaggedlist.extend(curxvals)
        flaggedlist.extend(xvals[:rmslen])
        flaggedlist.extend(xvals[-rmslen:])
        flaggedtimes = np.array(flaggedlist)
        notflagged = ~np.isin(xvals, flaggedtimes)

        raw_inroi = ((xvals - tpeak) > (-0.5 / 24)) & ((xvals - tpeak) < (1 / 24))
        baseline_mask = ((xvals - tpeak) > (-2 / 24)) & ((xvals - tpeak) < (-600 / 60 / 60 / 24)) & notflagged
        if not np.any(baseline_mask):
            baseline_mask = (~raw_inroi) & notflagged
        baseline = np.nanmedian(yvals[baseline_mask]) if np.any(baseline_mask) else 0.0
        yplot = yvals - baseline

        bbaseline_mask = ((bxvals - tpeak) > (-2 / 24)) & ((bxvals - tpeak) < (-600 / 60 / 60 / 24))
        if np.any(np.isfinite(byvals)):
            bbaseline = np.nanmedian(byvals[bbaseline_mask]) if np.any(bbaseline_mask) else np.nanmedian(byvals)
            byplot = byvals - bbaseline
        else:
            byplot = byvals

        fig = plt.figure(figsize=(10, 10))
        gs = GridSpec(6, 5, figure=fig)

        runrms = np.array(runrms)
        ax = fig.add_subplot(gs[-4, :])
        if runrms.size:
            ax.scatter(runrms[:, 0], runrms[:, 1])
        ax.axhline(args.rmsthresh)
        ax.set_yscale("log")

        ax = fig.add_subplot(gs[-3, :])
        ax.errorbar(xvals[notflagged], yplot[notflagged], yerr[notflagged], ls="none", color="black")
        ax.errorbar(xvals[~notflagged], yplot[~notflagged], yerr[~notflagged], ls="none", color="red", marker="X", label="flagged")
        ax.axvline(tpeak, ls="--", label="Peak")
        ax.set_ylabel("Difference flux (cts/s)")
        ax.legend()

        axb = fig.add_subplot(gs[-2, :])
        bnotflagged = ~np.isin(bxvals, flaggedtimes)
        axb.errorbar(bxvals[bnotflagged], byplot[bnotflagged], byerr[bnotflagged], ls="none", color="black")
        axb.errorbar(bxvals[~bnotflagged], byplot[~bnotflagged], byerr[~bnotflagged], ls="none", color="red", marker="X", label="flagged")
        axb.axvline(tpeak, ls="--", label="Peak")
        axb.set_ylabel("Bkg difference flux (cts/s)")
        axb.legend()

        inroi = ((xvals - tpeak) > (-0.5 / 24)) & ((xvals - tpeak) < (1 / 24))
        if not np.any(inroi):
            raise RuntimeError(f"No data near peak for {lcfile}")
        subxvals = xvals[inroi]
        subyvals = yplot[inroi]
        subxmag = subxvals[subyvals > 0]
        subymag = ctstomag(subyvals[subyvals > 0])
        otherflux = yplot[(~inroi) & notflagged]
        otherrms = np.sqrt(np.nanmean(otherflux**2)) if np.any(np.isfinite(otherflux)) else np.nan
        rms = np.sqrt(np.nanmean(subyvals[subyvals > 0] ** 2)) if np.any(subyvals > 0) else np.nan
        xmax = subxvals[np.nanargmin(np.abs(subxvals - tpeak))]
        ymax = subyvals[subxvals == xmax]
        errmax = yerr[inroi][subxvals == xmax]
        errmax = errmax[0] if errmax.size else np.nan
        ax.errorbar(xmax, ymax, errmax, ls="none", color="black", marker="*")

        peak_region = np.abs(xvals - tpeak) < (2 / 24)
        peak_flag_fraction = (
            np.sum(np.isin(xvals[peak_region], flaggedtimes)) / np.sum(peak_region)
            if np.sum(peak_region) > 0
            else 1.0
        )
        flagaroundmax = np.sum((flaggedtimes > (tpeak - (2 / 24))) & (flaggedtimes < (tpeak + (2 / 24))))
        prefix = ""
        if (
            (rms < (3 * otherrms))
            or (flagaroundmax > 0 and peak_flag_fraction > 0.5)
            or (otherrms > rms)
            or (np.abs(xmax - tpeak) * 24 * 60 * 60 > 400)
        ):
            prefix = "flag_"

        bkgsegment = ((xvals - tpeak) > (-2 / 24)) & ((xvals - tpeak) < (-600 / 60 / 60 / 24)) & notflagged
        noise_source = byplot if np.any(np.isfinite(byplot)) else yplot
        if np.any(bkgsegment):
            ctssensitivity = 3 * np.sqrt(np.nanmean(noise_source[bkgsegment] ** 2))
        elif np.any(np.isfinite(otherflux)):
            ctssensitivity = 3 * np.sqrt(np.nanmean(otherflux**2))
        else:
            ctssensitivity = np.nan
        magsensitivity = ctstomag(ctssensitivity)

        axmag = fig.add_subplot(gs[-1, :])
        axmag.invert_yaxis()
        axmag.axvline(tpeak, ls="--", label="Peak")
        maglimitx = xvals[yplot < ctssensitivity] if np.isfinite(ctssensitivity) else np.array([])
        maglimity = np.array([magsensitivity for _ in maglimitx])
        regmask = ((tpeak - subxmag) < (0.5 / 24)) & ((subxmag - tpeak) < 1 / 24)
        regxmag = subxmag[regmask]
        regymag = subymag[regmask]
        regflux = subyvals[subyvals > 0][regmask]
        detmask = regflux > (3 * otherrms) if np.isfinite(otherrms) else np.zeros_like(regymag, dtype=bool)
        detxmag = regxmag[detmask]
        detymag = regymag[detmask]
        axmag.scatter(xvals[(yplot > 0) & notflagged & (~np.isin(xvals, maglimitx))], ctstomag(yplot[(yplot > 0) & notflagged & (~np.isin(xvals, maglimitx))]), color="black")
        axmag.scatter(regxmag[~np.isin(regxmag, detxmag)], regymag[~np.isin(regxmag, detxmag)], color="blue", label="Near peak")
        axmag.scatter(detxmag, detymag, color="goldenrod", label=f"mag < {ctstomag(3 * otherrms):.2f}" if np.isfinite(otherrms) else "det")
        axmag.scatter(maglimitx, maglimity, color="black", marker="v", label=f"mag > {magsensitivity:.2f}" if np.isfinite(magsensitivity) else "limit")
        axmag.set_xlim(tpeak - (1 / 24), tpeak + (2 / 24))
        try:
            axmag.set_ylim(20, 0.9 * np.nanmin(detymag))
        except Exception:
            axmag.set_ylim(20, 10)
        axmag.set_xlabel("TJD (days)")
        axmag.set_ylabel("Difference Tmag")
        axmag.legend()

        numpts = len(detymag)
        if numpts == 0:
            prefix = "flag_"
        trigcand = Time(tpeak + 2457000, format="jd").iso
        title = "TESS Light Curve of " + lcfile + f" || RMS:{rms:9.2f} cts/s || {numpts} points || negatives: {np.sum(yplot[notflagged] < 0)}"

        imsize = 20
        detorbit_hint = None
        try:
            det_image, header, det_info = load_detection_images(tpeak, scale)
            detorbit_hint = str(det_info["detorbit"])
            fit_ra, fit_dec, fit_gal, ellip_ratio = fit_gaussian_panel(fig, gs, det_image, header, src, imsize)
            ra, dec, src_gal = fit_ra, fit_dec, fit_gal
            if (ellip_ratio > 1.5) or (ellip_ratio < 0.5):
                prefix = "flag_"
        except Exception:
            print(traceback.format_exc())

        add_reference_panel(fig, gs, scale, src, imsize)
        try:
            add_daily_rms_panels(fig, gs, candidate_day, scale, src, imsize, detorbit_hint=detorbit_hint)
        except Exception:
            print(traceback.format_exc())

        outfolder = os.path.join(os.path.dirname(os.path.dirname(lcfile)), "plots")
        os.makedirs(outfolder, exist_ok=True)
        outfile = os.path.basename(lcfile).replace(".txt", "") + f"cam{cam}_ccd{ccd}.jpg"
        if numpts == 1:
            prefix = prefix + "single_"
        fig.suptitle(
            f"{title}\nCLASS: {prefix} || (RA,Dec):({ra:7.4f},{dec:7.4f}) || "
            f"(col,row):({fcol:6.1f},{frow:6.1f})\nDatafile:{lcfile} || "
            f"(l,b): ({src_gal.l.deg:7.4f},{src_gal.b.deg:7.4f}) || peak: {trigcand}"
        )
        plt.savefig(os.path.join(outfolder, prefix + outfile))
        plt.close()

    except Exception:
        print(traceback.format_exc())


if __name__ == "__main__":
    import glob

    lcfiles = sorted(
        f
        for f in glob.glob(f"{args.photfile.replace('phot.data', '')}/lc/lc_outcatrms_*.txt")
        if "_cleaned" not in f and not f.endswith(".tmp")
    )
    if args.singleproc:
        for lc in tqdm(lcfiles):
            print(lc)
            plot(lc)
    else:
        with Pool(args.N) as pool:
            list(tqdm(pool.imap_unordered(plot, lcfiles), total=len(lcfiles)))
