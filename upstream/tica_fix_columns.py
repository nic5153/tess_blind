"""Apply the column cleaning (tica_common.clean_columns) to difference images
that were made before it existed, rewriting each slice's images.tar in place.
No downloads; conv_ members that already have COLCLEAN=T are left alone.

    python upstream/tica_fix_columns.py --sector 101 --cam 3 --ccd 3 [--orbit o1a ...]

After fixing an orbit, rebuild its daily RMS images: delete
$WORK_DIR/sectorN/camC_ccdD/oXX/ (and that cam/ccd's phot.data, lc/, bkg_phot/)
and rerun orbit_pipeline.sbatch.
"""
import argparse
import glob
import io
import os
import tarfile
import time
from multiprocessing import Pool

import numpy as np
from astropy.io import fits

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--orbit", nargs="*", default=None, help="default: every orbit present")
args = parser.parse_args()

ccd_dir = tc.ccd_dir(args.sector, args.cam, args.ccd)
orbits = args.orbit or sorted(os.path.basename(p) for p in glob.glob(f"{ccd_dir}/o[0-9][ab]"))


def fix_slice(slice_dir):
    tar_path = os.path.join(slice_dir, "images.tar")
    if not os.path.exists(tar_path):
        return slice_dir, 0, 0, "no images.tar"
    tmp = tar_path + ".fixing"
    nfix = nbad = 0
    with tarfile.open(tar_path, "r") as src, tarfile.open(tmp, "w") as dst:
        for member in src.getmembers():
            payload = src.extractfile(member).read()
            if member.name.startswith("conv_"):
                with fits.open(io.BytesIO(payload)) as hdul:
                    hdr = hdul[0].header.copy()
                    if not hdr.get("COLCLEAN", False):
                        diff = np.array(hdul[0].data, dtype=np.float32)
                        bad = tc.clean_columns(diff)
                        hdr["COLCLEAN"] = (True, "column medians removed, noisy columns zeroed")
                        hdr["NBADCOL"] = (len(bad), "columns zeroed by clean_columns")
                        hdr["BADCOLS"] = (tc.badcols_card(bad), "zeroed cols")
                        payload = tc.fits_bytes(diff, hdr)
                        nfix += 1
                        nbad += len(bad)
            info = tarfile.TarInfo(member.name)
            info.size = len(payload)
            info.mtime = int(time.time())
            dst.addfile(info, io.BytesIO(payload))
    if nfix:
        os.replace(tmp, tar_path)
    else:
        os.remove(tmp)
    return slice_dir, nfix, nbad, ""


for orbit in orbits:
    orbit_dir = os.path.join(ccd_dir, orbit)
    if os.path.exists(os.path.join(orbit_dir, ".purged")):
        print(f"{orbit}: purged, nothing to fix")
        continue
    slices = sorted(glob.glob(f"{orbit_dir}/slice*"))
    t0 = time.time()
    with Pool(tc.env_int("NPROC", 4)) as pool:
        results = pool.map(fix_slice, slices, chunksize=4)
    frames = sum(r[1] for r in results)
    cols = sum(r[2] for r in results)
    print(f"{orbit}: {frames} difference images cleaned in {len(slices)} slices "
          f"({cols / max(frames, 1):.1f} columns zeroed per image) in {time.time() - t0:.0f}s")
