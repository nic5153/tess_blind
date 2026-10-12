"""List the TICA FFIs for one sector/cam/ccd on MAST and assign each frame to
the HPCC-style half-orbit (o1a, o1b, o2a, o2b) and slice (sliceNNNN).

Half-orbits: each TICA orbit (o1/o2) is split at its largest internal cadence
gap (the mid-orbit downlink). Sectors without one are split at the midpoint
(HALF_SPLIT=midpoint) or kept whole as 'a' (HALF_SPLIT=none).

Slices: SLICE_SIZE frames per slice, or SLICE_SIZE=auto to reproduce the HPCC
layout, which used round(N_halforbit / 450) frames per slice.

Writes $DATA_DIR/s00XX/camC-ccdD/manifest.tsv. Downloads nothing but the listing.
"""
import argparse
import os
import urllib.request
from collections import defaultdict

import numpy as np

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--refresh", action="store_true", help="rebuild an existing manifest")
args = parser.parse_args()

outfile = tc.manifest_path(args.sector, args.cam, args.ccd)
if os.path.exists(outfile) and not args.refresh:
    print(f"Manifest exists: {outfile} (use --refresh to rebuild)")
    raise SystemExit(0)

slice_size = tc.env_str("SLICE_SIZE", "auto")
slice_target = tc.env_int("SLICE_TARGET", 450)
gap_min = tc.env_int("HALF_GAP_MIN_CADENCES", 20)
half_split = tc.env_str("HALF_SPLIT", "midpoint")

url = f"{tc.TICA_URL}/s{args.sector:04d}/cam{args.cam}-ccd{args.ccd}/"
print(f"Listing {url}")
with urllib.request.urlopen(url, timeout=120) as resp:
    html = resp.read().decode("utf-8", errors="replace")

frames = {}
for m in tc.NAME_RE.finditer(html):
    sec, torbit, cad, cam, ccd = m.groups()
    if int(sec) == args.sector and int(cam) == args.cam and int(ccd) == args.ccd:
        frames[m.group(0)] = (torbit, int(cad))

if not frames:
    raise SystemExit(f"ERROR: no TICA frames found at {url}")

by_orbit = defaultdict(list)
for name, (torbit, cad) in frames.items():
    by_orbit[torbit].append((cad, name))

rows = []
for torbit in sorted(by_orbit):
    entries = sorted(by_orbit[torbit])
    cads = np.array([c for c, _ in entries])
    gaps = np.diff(cads) - 1
    split = None
    if gaps.size and gaps.max() >= gap_min:
        split = int(np.argmax(gaps)) + 1
        how = f"gap of {gaps.max()} cadences"
    elif half_split == "midpoint" and len(entries) > 1:
        split = len(entries) // 2
        how = "midpoint (no mid-orbit gap)"
    else:
        how = "not split"

    halves = {"a": entries[:split] if split else entries, "b": entries[split:] if split else []}
    for half, members in halves.items():
        if not members:
            continue
        orbit = f"{torbit}{half}"
        n = len(members)
        per = max(1, round(n / slice_target)) if slice_size == "auto" else max(1, int(slice_size))
        nslices = (n + per - 1) // per
        print(f"  {torbit} -> {orbit}: {n} frames, {per} per slice, {nslices} slices [{how}]")
        for i, (cad, name) in enumerate(members):
            rows.append((name, torbit, cad, orbit, f"slice{i // per:04d}"))

os.makedirs(os.path.dirname(outfile), exist_ok=True)
tmp = outfile + ".tmp"
with open(tmp, "w") as f:
    f.write("name\ttica_orbit\tcadence\torbit\tslice\n")
    for row in rows:
        f.write("\t".join(str(v) for v in row) + "\n")
os.replace(tmp, outfile)
print(f"Wrote {len(rows)} frames to {outfile}")
