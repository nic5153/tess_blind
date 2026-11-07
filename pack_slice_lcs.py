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
PARTIAL = "lc_partial.tsv"  # incomplete slice: kept, but re-measured on resume
_expected = {}


def expected_for(slice_dir):
    """Candidate names (from the orbit's phot.data) and frame count (dates) the
    slice's photometry must contain to count as complete."""
    orbit_dir = os.path.dirname(slice_dir)
    if orbit_dir not in _expected:
        names = set()
        pf = os.path.join(orbit_dir, "phot.data")
        if os.path.exists(pf):
            with open(pf) as f:
                for line in f:
                    parts = line.split()
                    if len(parts) == 6:
                        names.add(os.path.basename(parts[4]))
        _expected[orbit_dir] = names
    with open(os.path.join(slice_dir, "dates")) as f:
        nframes = sum(1 for line in f if line.strip())
    return _expected[orbit_dir], nframes


def pack(d):
    pieces = sorted(f for f in os.listdir(d)
                    if f.startswith("lc_") and f not in (COMBINED, PARTIAL) and not f.endswith(".part"))
    if not pieces:
        return d, 0, 0, "ok"
    slice_dir = d[: -len(os.path.join("bkg_phot", args.lcdir))].rstrip(os.sep) \
        if d.endswith(os.path.join("bkg_phot", args.lcdir)) else os.path.dirname(d)
    want_names, nframes = expected_for(slice_dir)
    complete = bool(want_names) and set(pieces) == want_names
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
        rows = [row for row in text.splitlines() if row.strip()]
        complete = complete and len(rows) == nframes
        lines.extend(f"{name} {row}\n" for row in rows)
    status = "packed" if complete else "partial"
    if args.dry_run:
        return d, len(pieces), len(lines), f"dry-run-{status}"

    # Reuse an existing inode: the combined file if present, else the first piece.
    # An incomplete slice goes to lc_partial.tsv, which collect_phot.py ignores
    # and a PHOT_RESUME run re-measures.
    target_path = combined_path if complete else os.path.join(d, PARTIAL)
    if complete and os.path.exists(combined_path):
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
        os.replace(holder, target_path)
    for name in pieces:
        if name != renamed:
            os.remove(os.path.join(d, name))
    return d, len(pieces), len(lines), status


pattern = os.path.join(args.data_dir, f"s{args.sector:04d}", "cam*-ccd*", "o[0-9][ab]", "slice[0-9]*")
dirs = []
for sdir in sorted(glob.glob(pattern)):
    for sub in (args.lcdir, os.path.join("bkg_phot", args.lcdir)):
        if os.path.isdir(os.path.join(sdir, sub)):
            dirs.append(os.path.join(sdir, sub))
print(f"{len(dirs)} slice light-curve directories under {args.data_dir}/s{args.sector:04d}"
      + (" (dry run)" if args.dry_run else ""))

pieces = rows = failed = partial = 0
with Pool(args.workers) as pool:
    for i, (d, npieces, nrows, status) in enumerate(pool.imap_unordered(pack, dirs, chunksize=8), 1):
        pieces += npieces
        rows += nrows
        if status.startswith("VERIFY"):
            failed += 1
            print(f"  {d}: {status}")
        elif status.endswith("partial"):
            partial += 1
        if i % 2000 == 0:
            print(f"  {i}/{len(dirs)} directories, {pieces} pieces so far", flush=True)
verb = "would pack" if args.dry_run else "packed"
print(f"{verb} {pieces} piece files ({rows} rows) into at most {len(dirs)} files; "
      f"{partial} directories incomplete (kept as {PARTIAL}, re-measured on resume); "
      f"{failed} failed verification and were left untouched")
