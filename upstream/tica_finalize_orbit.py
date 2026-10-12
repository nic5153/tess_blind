"""Write $DATA_DIR/s00XX/camC-ccdD/oXX/rms_list_tica: the frames make_rms.py may
bin into daily RMS images, one name per line in time order.

A frame is left out if
  * its QUAL_BIT has any RMS_QUAL_MASK bit set (default -1: any TICA flag), or
  * its difference image is an outlier (bad frames TICA does not flag, e.g. a
    pointing excursion; one such frame in a day's RMS image yields thousands of
    false single-day candidates):
      - residual scatter conv_std > RMS_MAX_NOISE_RATIO (default 3) x the orbit
        median, or
      - kernel flux scale |kscale - 1| > RMS_MAX_KSCALE_DEV (default 0.1).

Excluded frames and reasons go to oXX/rms_excluded.txt. Frames stay in the
slice tarballs and light curves; only the RMS images skip them.
"""
import argparse
import glob
import os

import numpy as np

import tica_common as tc

parser = argparse.ArgumentParser()
parser.add_argument("--sector", type=int, required=True)
parser.add_argument("--cam", type=int, required=True)
parser.add_argument("--ccd", type=int, required=True)
parser.add_argument("--orbit", required=True)
args = parser.parse_args()

mask = tc.env_int("RMS_QUAL_MASK", -1)
max_noise = tc.env_float("RMS_MAX_NOISE_RATIO", 3.0)
max_kdev = tc.env_float("RMS_MAX_KSCALE_DEV", 0.1)
orbit_dir = os.path.join(tc.ccd_dir(args.sector, args.cam, args.ccd), args.orbit)

expected = {r["slice"] for r in tc.read_manifest(args.sector, args.cam, args.ccd)
            if r["orbit"] == args.orbit}
done = {os.path.basename(os.path.dirname(p)) for p in glob.glob(f"{orbit_dir}/slice*/.done")}
missing = sorted(expected - done)

frames = {}  # name -> [tjd, qual, conv_std, kscale]
for slice_name in sorted(done):
    sdir = os.path.join(orbit_dir, slice_name)
    with open(os.path.join(sdir, ".quality")) as f:
        for line in f:
            name, qual, tjd = line.split()
            frames[name] = [float(tjd), int(qual), np.nan, np.nan]
    stats = os.path.join(sdir, ".stats.tsv")
    if os.path.exists(stats):
        with open(stats) as f:
            header = f.readline().rstrip("\n").split("\t")
            for line in f:
                row = dict(zip(header, line.rstrip("\n").split("\t")))
                if row.get("name") in frames:
                    frames[row["name"]][2] = float(row["conv_std"])
                    frames[row["name"]][3] = float(row["kscale"])

noise = np.array([v[2] for v in frames.values()])
median_noise = float(np.nanmedian(noise)) if np.isfinite(noise).any() else np.nan

kept, excluded = [], []
for name, (tjd, qual, conv_std, kscale) in frames.items():
    reasons = []
    if qual & mask:
        reasons.append(f"QUAL_BIT={qual}")
    if np.isfinite(conv_std) and conv_std > max_noise * median_noise:
        reasons.append(f"conv_std={conv_std:.0f}>{max_noise:g}x{median_noise:.0f}")
    if np.isfinite(kscale) and abs(kscale - 1.0) > max_kdev:
        reasons.append(f"kscale={kscale:.3f}")
    (excluded if reasons else kept).append((tjd, name, "; ".join(reasons)))
kept.sort()
excluded.sort()

out = os.path.join(orbit_dir, "rms_list_tica")
with open(out + ".tmp", "w") as f:
    for _, name, _ in kept:
        f.write(name + "\n")
os.replace(out + ".tmp", out)
with open(os.path.join(orbit_dir, "rms_excluded.txt"), "w") as f:
    for tjd, name, why in excluded:
        f.write(f"{name} {tjd:.6f} {why}\n")

print(f"{out}: {len(kept)} frames listed, {len(excluded)} excluded "
      f"(median conv_std {median_noise:.0f}; see rms_excluded.txt)")
for tjd, name, why in excluded[:10]:
    print(f"  excluded TJD {tjd:.4f} {name}: {why}")
if missing:
    print(f"WARNING: {len(missing)} of {len(expected)} slices not processed: "
          f"{', '.join(missing[:10])}{' ...' if len(missing) > 10 else ''}")
    raise SystemExit(1)
