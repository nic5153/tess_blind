"""Gather the per-slice photometry written by phot_slice.py into the usual
per-candidate light curves (the job copy_phot_em2.sbatch does with cat >> for
ISIS phot2 output):

    <orbit>/sliceNNNN/<lcdir>/lc_combined.tsv
        -> <outdir>/<lcdir>/lc_<candidate>.txt
    <orbit>/sliceNNNN/bkg_phot/<lcdir>/lc_combined.tsv
        -> <outdir>/bkg_phot/<lcdir>/lc_<candidate>.txt

Rows are appended (light curves accumulate over orbits, as before);
orbit_pipeline.sbatch then sorts and de-duplicates them by time.

    python collect_phot.py --orbit-dir $DATA_DIR/s0102/cam3-ccd3/o1a \
        --outdir sector102/cam3_ccd3 [--lcdir lc]
"""
import argparse
import glob
import os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("--orbit-dir", required=True)
parser.add_argument("--outdir", required=True, help="sector*/camC_ccdD directory")
parser.add_argument("--lcdir", default="lc")
parser.add_argument("--flush-every", type=int, default=50, help="slices buffered between writes")
args = parser.parse_args()

slices = sorted(glob.glob(os.path.join(args.orbit_dir, "slice[0-9]*")))
for sub in (args.lcdir, os.path.join("bkg_phot", args.lcdir)):
    target = os.path.join(args.outdir, sub)
    os.makedirs(target, exist_ok=True)
    buffer = defaultdict(list)
    nfiles = nrows = missing = 0

    def flush():
        global nrows
        for name, lines in buffer.items():
            with open(os.path.join(target, name), "a") as out:
                out.writelines(lines)
            nrows += len(lines)
        buffer.clear()

    for i, sdir in enumerate(slices, 1):
        path = os.path.join(sdir, sub, "lc_combined.tsv")
        if not os.path.exists(path):
            missing += 1
            continue
        nfiles += 1
        with open(path) as f:
            for line in f:
                name, _, rest = line.partition(" ")
                if rest:
                    buffer[name].append(rest)
        if i % args.flush_every == 0:
            flush()
    flush()
    print(f"collect_phot: {sub}: {nfiles} slices, {nrows} rows -> {target}"
          + (f" ({missing} slices had no photometry)" if missing else ""))
