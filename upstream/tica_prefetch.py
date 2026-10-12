"""Download TICA frames ahead of time into $RAW_CACHE/s00XX/camC-ccdD/pool/, for
clusters whose compute nodes cannot reach the internet. Run it on a login or
data-transfer node; the processing jobs then take frames from the pool (and
delete them as usual) instead of downloading.

    python upstream/tica_prefetch.py --sector 101 --cam 3 --ccd 3 --orbit o2a [--threads 8]

Frames of slices that are already done are skipped. If ref.fits does not exist
yet, the reference candidates (from REF_ORBIT, default the first half-orbit)
are fetched as well. Each CCD half-orbit is ~15-50 GB of raw FFIs.
"""
import argparse
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import numpy as np

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--orbit", nargs="+", required=True)
parser.add_argument("--threads", type=int, default=8)
args = parser.parse_args()

here = os.path.dirname(os.path.abspath(__file__))
if not os.path.exists(tc.manifest_path(args.sector, args.cam, args.ccd)):
    subprocess.run([sys.executable, os.path.join(here, "tica_manifest.py"), "--sector", str(args.sector),
                    "--cam", str(args.cam), "--ccd", str(args.ccd)], check=True)
rows = sorted(tc.read_manifest(args.sector, args.cam, args.ccd), key=lambda r: r["cadence"])
ccd_dir = tc.ccd_dir(args.sector, args.cam, args.ccd)

wanted = []
for r in rows:
    if r["orbit"] in args.orbit and not os.path.exists(os.path.join(ccd_dir, r["orbit"], r["slice"], ".done")):
        wanted.append(r["name"])

# Reference candidates, chosen exactly as tica_make_ref.py chooses them.
if not os.path.exists(os.path.join(ccd_dir, "ref.fits")):
    ref_orbit = os.environ.get("REF_ORBIT") or sorted({r["orbit"] for r in rows})[0]
    ref_rows = [r for r in rows if r["orbit"] == ref_orbit]
    ncand = max(tc.env_int("REF_NFRAMES", 20), tc.env_int("REF_NCAND", 60))
    picks = np.unique(np.linspace(0, len(ref_rows) - 1, min(ncand, len(ref_rows))).astype(int))
    wanted += [ref_rows[i]["name"] for i in picks]

pool = tc.pool_dir(args.sector, args.cam, args.ccd)
os.makedirs(pool, exist_ok=True)
todo = sorted({n for n in wanted if not os.path.exists(os.path.join(pool, n))})
print(f"Prefetching {len(todo)} frames (~{len(todo) * 17.8 / 1024:.0f} GB) into {pool} "
      f"with {args.threads} threads")


def get(name):
    tc.download(tc.frame_url(args.sector, args.cam, args.ccd, name), os.path.join(pool, name))
    return name


done = 0
with ThreadPoolExecutor(args.threads) as ex:
    for _ in ex.map(get, todo):
        done += 1
        if done % 200 == 0 or done == len(todo):
            print(f"  {done}/{len(todo)}", flush=True)
print("Prefetch complete.")
