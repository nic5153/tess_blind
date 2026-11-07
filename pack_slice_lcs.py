"""Losslessly repack per-slice light-curve pieces into one file per slice.

The first phot_slice.py wrote one small file per candidate per slice
(sliceNNNN/lc/lc_*.txt and sliceNNNN/bkg_phot/lc/lc_*.txt), which can exceed an
HPCC file-count quota. This packs each such directory into the current format,
a single lc_combined.tsv ("<lc file name> <row>" per line), that
collect_phot.py merges into the per-candidate light curves.

No data is dropped and no new file is needed (it works over an inode quota):
the combined content is written into one of the existing piece files, which is
then renamed; the row count is verified before any other piece is deleted, and
on a mismatch the directory is restored and left untouched.

    python pack_slice_lcs.py --sector 102 --dry-run      # report only
    python pack_slice_lcs.py --sector 102 --workers 8    # pack
"""
import argparse
import glob
import os
from multiprocessing import Pool

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--data-dir", default=os.environ.get("DATA_DIR", ""))
parser.add_argument("--lcdir", default="lc")
parser.add_argument("--workers", type=int, default=4)
parser.add_argument("--dry-run", action="store_true")
args = parser.parse_args()
if not args.data_dir:
    raise SystemExit("ERROR: pass --data-dir or source setup_vars_hpcc.sh (DATA_DIR)")

COMBINED = "lc_combined.tsv"


def pack(d):
    pieces = sorted(f for f in os.listdir(d)
                    if f.startswith("lc_") and f != COMBINED and not f.endswith(".part"))
    if not pieces:
        return d, 0, 0, "ok"
    combined_path = os.path.join(d, COMBINED)
    lines = []
    if os.path.exists(combined_path):
        with open(combined_path) as f:
            lines.extend(f.readlines())
    originals = {}
    for name in pieces:
        with open(os.path.join(d, name)) as f:
            text = f.read()
        originals[name] = text
        for row in text.splitlines():
            if row.strip():
                lines.append(f"{name} {row}\n")
    if args.dry_run:
        return d, len(pieces), len(lines), "dry-run"

    # Reuse an existing inode: the combined file if present, else the first piece.
    if os.path.exists(combined_path):
        holder, renamed = combined_path, None
    else:
        holder, renamed = os.path.join(d, pieces[0]), pieces[0]
    with open(holder, "w") as f:
        f.writelines(lines)
    with open(holder) as f:
        written = sum(1 for _ in f)
    if written != len(lines):
        if renamed is not None:  # restore the reused piece
            with open(holder, "w") as f:
                f.write(originals[renamed])
        return d, len(pieces), len(lines), f"VERIFY FAILED ({written} != {len(lines)}); left as is"
    if renamed is not None:
        os.replace(holder, combined_path)
    for name in pieces:
        if name != renamed:
            os.remove(os.path.join(d, name))
    return d, len(pieces), len(lines), "packed"


pattern = os.path.join(args.data_dir, f"s{args.sector:04d}", "cam*-ccd*", "o[0-9][ab]", "slice[0-9]*")
dirs = []
for sdir in sorted(glob.glob(pattern)):
    for sub in (args.lcdir, os.path.join("bkg_phot", args.lcdir)):
        if os.path.isdir(os.path.join(sdir, sub)):
            dirs.append(os.path.join(sdir, sub))
print(f"{len(dirs)} slice light-curve directories under {args.data_dir}/s{args.sector:04d}"
      + (" (dry run)" if args.dry_run else ""))

pieces = rows = failed = 0
with Pool(args.workers) as pool:
    for i, (d, npieces, nrows, status) in enumerate(pool.imap_unordered(pack, dirs, chunksize=8), 1):
        pieces += npieces
        rows += nrows
        if status.startswith("VERIFY"):
            failed += 1
            print(f"  {d}: {status}")
        if i % 2000 == 0:
            print(f"  {i}/{len(dirs)} directories, {pieces} pieces so far", flush=True)
verb = "would pack" if args.dry_run else "packed"
print(f"{verb} {pieces} piece files ({rows} rows) into at most {len(dirs)} combined files; "
      f"{failed} directories failed verification and were left untouched")
