#!/usr/bin/env python

import argparse
import csv
import glob
import os
import re
from pathlib import Path

import numpy as np
from tess_stars2px import tess_stars2px_reverse_function_entry


def parse_args():
    parser = argparse.ArgumentParser(
        description="Export sector*/multiple plot candidates to a CSV for source follow-up."
    )
    parser.add_argument(
        "--root",
        default=".",
        help="Pipeline root containing sector*/ directories. Default: current directory.",
    )
    parser.add_argument(
        "--outfile",
        default="multiple_candidates.csv",
        help="Output CSV path. Default: multiple_candidates.csv",
    )
    parser.add_argument(
        "--html",
        default="multiple_candidates.html",
        help="Output HTML review page. Default: multiple_candidates.html",
    )
    return parser.parse_args()


def parse_sector(path):
    match = re.search(r"sector(\d+)", path.replace("\\", "/"))
    if match is None:
        return None
    return int(match.group(1))


def parse_camccd(path):
    match = re.search(r"cam(\d+)_ccd(\d+)", path)
    if match is None:
        return None, None
    return int(match.group(1)), int(match.group(2))


def parse_rms_day(path):
    match = re.search(r"outcatrms_day_(\d+)", os.path.basename(path))
    if match is None:
        return ""
    return match.group(1)


def source_lc_candidates_from_plot(plot_path, cam, ccd):
    stem = Path(plot_path).stem
    suffix = f"cam{cam}_ccd{ccd}"
    if stem.endswith(suffix):
        stem = stem[: -len(suffix)]
    for prefix in ("flag_single_", "flag_", "single_"):
        if stem.startswith(prefix):
            stem = stem[len(prefix) :]

    raw_stem = stem[:-8] if stem.endswith("_cleaned") else stem
    stems = [stem]
    if raw_stem != stem:
        stems.append(raw_stem)

    candidates = []
    for candidate_stem in stems:
        candidates.extend(
            [
                f"lc/{candidate_stem}.txt",
                f"lcGRB/{candidate_stem}.txt",
                f"{candidate_stem}.txt",
                candidate_stem,
            ]
        )
    return candidates


def load_phot_data(photfile):
    if not os.path.exists(photfile) or os.path.getsize(photfile) == 0:
        return None
    data = np.loadtxt(
        photfile,
        dtype=[
            ("fcol", "f8"),
            ("frow", "f8"),
            ("icol", "i8"),
            ("irow", "i8"),
            ("fname", "<U256"),
            ("flag", "i4"),
        ],
    )
    if data.size == 1:
        data = data[np.newaxis]
    return data


def find_phot_row(phot_data, lc_candidates):
    if phot_data is None:
        return None
    candidates = set(lc_candidates)
    candidates.update(os.path.basename(name) for name in lc_candidates)
    for row in phot_data:
        if str(row["fname"]) in candidates:
            return row
    return None


def find_orbit(sector_dir, cam, ccd, rms_day):
    if not rms_day:
        return ""
    pattern = os.path.join(sector_dir, f"cam{cam}_ccd{ccd}", "o??", f"rms_day_{rms_day}.fits")
    matches = sorted(glob.glob(pattern))
    if not matches:
        return ""
    return os.path.basename(os.path.dirname(matches[0]))


def main():
    args = parse_args()
    root = os.path.abspath(args.root)
    plot_files = sorted(glob.glob(os.path.join(root, "sector*", "multiple", "*.jpg")))

    rows = []
    scinfo_cache = {}
    phot_cache = {}

    for plot_file in plot_files:
        sector = parse_sector(plot_file)
        cam, ccd = parse_camccd(os.path.basename(plot_file))
        if sector is None or cam is None or ccd is None:
            print(f"WARNING: could not parse sector/cam/ccd from {plot_file}")
            continue

        sector_dir = os.path.join(root, f"sector{sector}")
        photfile = os.path.join(sector_dir, f"cam{cam}_ccd{ccd}", "phot.data")
        if photfile not in phot_cache:
            phot_cache[photfile] = load_phot_data(photfile)

        lc_candidates = source_lc_candidates_from_plot(plot_file, cam, ccd)
        phot_row = find_phot_row(phot_cache[photfile], lc_candidates)
        if phot_row is None:
            print(f"WARNING: no phot.data match for {plot_file} as {lc_candidates[0]}")
            fcol = frow = ra = dec = ""
            lc_file = lc_candidates[0]
        else:
            fcol = float(phot_row["fcol"])
            frow = float(phot_row["frow"])
            lc_file = str(phot_row["fname"])
            key = (sector, cam, ccd)
            scinfo = scinfo_cache.get(key)
            ra, dec, scinfo = tess_stars2px_reverse_function_entry(
                sector, cam, ccd, fcol, frow, scInfo=scinfo
            )
            scinfo_cache[key] = scinfo

        rms_day = parse_rms_day(plot_file)
        orbit = find_orbit(sector_dir, cam, ccd, rms_day)

        rows.append(
            {
                "sector": sector,
                "cam_ccd": f"cam{cam}_ccd{ccd}",
                "orbit": orbit,
                "rms_date": rms_day,
                "ra": ra,
                "dec": dec,
                "suspected source type": "",
                "notes": "",
                "plot_file": os.path.relpath(plot_file, root).replace("\\", "/"),
                "lc_file": lc_file,
                "fcol": fcol,
                "frow": frow,
            }
        )

    fieldnames = [
        "sector",
        "cam_ccd",
        "orbit",
        "rms_date",
        "ra",
        "dec",
        "suspected source type",
        "notes",
        "plot_file",
        "lc_file",
        "fcol",
        "frow",
    ]
    with open(args.outfile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    write_html(args.html, rows)
    print(f"Wrote {len(rows)} candidates to {args.outfile}")
    print(f"Wrote HTML review page to {args.html}")


def html_escape(value):
    text = "" if value is None else str(value)
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def write_html(outfile, rows):
    rows = sorted(rows, key=lambda r: (int(r["sector"]), r["cam_ccd"], r["orbit"], r["rms_date"], r["plot_file"]))
    sectors = []
    for row in rows:
        if row["sector"] not in sectors:
            sectors.append(row["sector"])

    html = """<html>
<! -- -->

<head>
<title>Multiple Candidates</title>
<style>
body { background: #ffffff; font-family: sans-serif; }
img { max-width: 1000px; width: 100%; height: auto; }
table { border-collapse: collapse; margin-bottom: 1em; }
td, th { border: 1px solid #cccccc; padding: 4px 6px; font-size: 13px; }
figure { margin-bottom: 2em; }
figcaption { font-size: 14px; }
</style>
</head>
<body bgcolor="#ffffff">
<h2>Multiple Candidates</h2>
"""

    for sector in sectors:
        sector_rows = [r for r in rows if r["sector"] == sector]
        html += f"<h3>Sector {html_escape(sector)} ({len(sector_rows)} candidates)</h3>\n"
        for row in sector_rows:
            html += "<figure>\n"
            html += f"""<figcaption>
<table>
<tr><th>sector</th><th>cam_ccd</th><th>orbit</th><th>rms_date</th><th>ra</th><th>dec</th><th>suspected source type</th><th>notes</th></tr>
<tr>
<td>{html_escape(row["sector"])}</td>
<td>{html_escape(row["cam_ccd"])}</td>
<td>{html_escape(row["orbit"])}</td>
<td>{html_escape(row["rms_date"])}</td>
<td>{html_escape(row["ra"])}</td>
<td>{html_escape(row["dec"])}</td>
<td>{html_escape(row["suspected source type"])}</td>
<td>{html_escape(row["notes"])}</td>
</tr>
</table>
{html_escape(row["plot_file"])}
</figcaption>\n"""
            html += f"""  <img src="{html_escape(row["plot_file"])}">\n"""
            html += "</figure>\n"

    html += "</body>\n</html>\n"
    with open(outfile, "w", newline="") as f:
        f.write(html)


if __name__ == "__main__":
    main()
