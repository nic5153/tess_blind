#!/usr/bin/env python

import argparse
import csv
import glob
import os
import re
from pathlib import Path

import numpy as np

try:
    from astropy.coordinates import SkyCoord
    import astropy.units as astrou
except ImportError:
    SkyCoord = None
    astrou = None

try:
    from tess_stars2px import tess_stars2px_reverse_function_entry
except ImportError:
    tess_stars2px_reverse_function_entry = None


BASE_FIELDNAMES = [
    "sector",
    "cam_ccd",
    "orbit",
    "rms_date",
    "ra",
    "dec",
    "l",
    "b",
    "suspected source type",
    "notes",
    "plot_file",
    "lc_file",
    "fcol",
    "frow",
]

CATALOG_FIELDNAMES = [
    "catalog_status",
    "catalog_matches",
    "color_summary",
    "gaia_source_id",
    "gaia_dist_arcsec",
    "gaia_parallax",
    "gaia_g_mag",
    "gaia_bp_mag",
    "gaia_rp_mag",
    "gaia_bp_rp",
    "panstarrs_dist_arcsec",
    "panstarrs_g_mag",
    "panstarrs_r_mag",
    "panstarrs_i_mag",
    "panstarrs_z_mag",
    "panstarrs_y_mag",
    "panstarrs_g_r",
    "panstarrs_r_i",
    "des_dist_arcsec",
    "des_g_mag",
    "des_r_mag",
    "des_i_mag",
    "des_z_mag",
    "des_y_mag",
    "des_g_r",
    "des_r_i",
    "decaps_match_count",
    "decaps_dist_arcsec",
    "decaps_g_mag",
    "decaps_r_mag",
    "decaps_i_mag",
    "decaps_z_mag",
    "decaps_y_mag",
    "decaps_g_r",
    "decaps_r_i",
]

EDITABLE_FIELDNAMES = ["suspected source type", "notes"]


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
    parser.add_argument(
        "--enriched-csv",
        default="multiple_candidates_enriched.csv",
        help=(
            "Optional catalog-enriched CSV to merge into the CSV/HTML output "
            "when it exists. Default: multiple_candidates_enriched.csv"
        ),
    )
    parser.add_argument(
        "--plot-root",
        action="append",
        default=[],
        help=(
            "Additional root containing a multiple/ plot directory. "
            "Can be used more than once."
        ),
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


def collect_plot_files(root, extra_roots):
    plot_files = []
    plot_files.extend(glob.glob(os.path.join(root, "sector*", "multiple", "*.jpg")))

    for plot_root in extra_roots:
        abs_plot_root = os.path.abspath(plot_root)
        plot_files.extend(glob.glob(os.path.join(abs_plot_root, "multiple", "*.jpg")))

    return sorted(dict.fromkeys(plot_files))


def html_image_path(plot_file, html_file):
    html_dir = os.path.dirname(os.path.abspath(html_file)) or os.getcwd()
    return os.path.relpath(plot_file, html_dir).replace("\\", "/")


def radec_to_galactic(ra, dec):
    if SkyCoord is None or astrou is None:
        return "", ""
    try:
        if ra == "" or dec == "":
            return "", ""
        sc = SkyCoord(float(ra), float(dec), unit="deg")
    except (TypeError, ValueError):
        return "", ""
    return float(sc.galactic.l.deg), float(sc.galactic.b.deg)


def blank_catalog_fields():
    return {field: "" for field in CATALOG_FIELDNAMES}


def row_key(row):
    plot_file = str(row.get("plot_file", "")).replace("\\", "/")
    if plot_file:
        return ("plot_file", plot_file)
    return (
        "coords",
        str(row.get("sector", "")),
        str(row.get("cam_ccd", "")),
        str(row.get("orbit", "")),
        str(row.get("rms_date", "")),
        str(row.get("ra", "")),
        str(row.get("dec", "")),
    )


def load_existing_metadata(paths):
    metadata = {}
    for path in paths:
        if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                metadata[row_key(row)] = row
    return metadata


def merge_existing_metadata(rows, paths):
    metadata = load_existing_metadata(paths)
    preserve_fields = EDITABLE_FIELDNAMES + CATALOG_FIELDNAMES
    for row in rows:
        previous = metadata.get(row_key(row))
        if previous is None:
            continue
        for field in preserve_fields:
            if field in previous and previous[field] != "":
                row[field] = previous[field]


def main():
    args = parse_args()
    root = os.path.abspath(args.root)
    plot_files = collect_plot_files(root, args.plot_root)

    rows = []
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
            if tess_stars2px_reverse_function_entry is None:
                ra = dec = ""
            else:
                # Match complicatedplot.py exactly so CSV/HTML coordinates
                # agree with the RA/Dec printed on the plot image.
                ra, dec, _ = tess_stars2px_reverse_function_entry(sector, cam, ccd, fcol, frow)

        gal_l, gal_b = radec_to_galactic(ra, dec)
        rms_day = parse_rms_day(plot_file)
        orbit = find_orbit(sector_dir, cam, ccd, rms_day)

        row = {
            "sector": sector,
            "cam_ccd": f"cam{cam}_ccd{ccd}",
            "orbit": orbit,
            "rms_date": rms_day,
            "ra": ra,
            "dec": dec,
            "l": gal_l,
            "b": gal_b,
            "suspected source type": "",
            "notes": "",
            "plot_file": html_image_path(plot_file, args.html),
            "lc_file": lc_file,
            "fcol": fcol,
            "frow": frow,
        }
        row.update(blank_catalog_fields())
        rows.append(row)

    merge_existing_metadata(rows, [args.outfile, args.enriched_csv])

    fieldnames = BASE_FIELDNAMES + CATALOG_FIELDNAMES
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


def html_catalog_summary(row):
    parts = []
    if row.get("catalog_matches"):
        parts.append(row["catalog_matches"])
    if row.get("catalog_status") and row["catalog_status"] != "ok":
        parts.append(row["catalog_status"])
    return "; ".join(parts)


def html_color_summary(row):
    if row.get("color_summary"):
        return row["color_summary"]
    parts = []
    if row.get("gaia_bp_rp"):
        parts.append(f"Gaia BP-RP={row['gaia_bp_rp']}")
    if row.get("panstarrs_g_r"):
        parts.append(f"PS1 g-r={row['panstarrs_g_r']}")
    if row.get("des_g_r"):
        parts.append(f"DES g-r={row['des_g_r']}")
    if row.get("decaps_g_r"):
        parts.append(f"DECaPS g-r={row['decaps_g_r']}")
    return "; ".join(parts)


def write_html(outfile, rows):
    rows = sorted(rows, key=lambda r: (int(r["sector"]), r["cam_ccd"], r["orbit"], r["rms_date"], r["plot_file"]))
    sectors = []
    for row in rows:
        if row["sector"] not in sectors:
            sectors.append(row["sector"])

    total = len(rows)
    nav_links = " ".join(
        f"""<a href="#sector-{html_escape(sector)}">Sector {html_escape(sector)}</a>"""
        for sector in sectors
    )

    html = """<html>
<! -- -->

<head>
<title>Multiple Candidates</title>
<style>
body { background: #ffffff; font-family: sans-serif; }
nav { position: sticky; top: 0; background: #ffffff; border-bottom: 1px solid #cccccc; padding: 8px 0; margin-bottom: 1em; }
nav a { display: inline-block; margin-right: 10px; padding: 4px 8px; border: 1px solid #cccccc; color: #000000; text-decoration: none; font-size: 14px; }
img { max-width: 1000px; width: 100%; height: auto; }
table { border-collapse: collapse; margin-bottom: 1em; }
td, th { border: 1px solid #cccccc; padding: 4px 6px; font-size: 13px; }
section { break-before: page; page-break-before: always; margin-bottom: 3em; }
section:first-of-type { break-before: auto; page-break-before: auto; }
figure { margin-bottom: 2em; }
figcaption { font-size: 14px; }
</style>
</head>
<body bgcolor="#ffffff">
"""
    html += f"<h2>Multiple Candidates ({total})</h2>\n"
    if nav_links:
        html += f"<nav>{nav_links}</nav>\n"

    for sector in sectors:
        sector_rows = [r for r in rows if r["sector"] == sector]
        html += f"""<section id="sector-{html_escape(sector)}">
<h3>Sector {html_escape(sector)} ({len(sector_rows)} candidates)</h3>
"""
        for row in sector_rows:
            html += "<figure>\n"
            html += f"""<figcaption>
<table>
<tr><th>sector</th><th>cam_ccd</th><th>orbit</th><th>rms_date</th><th>ra</th><th>dec</th><th>l</th><th>b</th><th>suspected source type</th><th>notes</th></tr>
<tr>
<td>{html_escape(row["sector"])}</td>
<td>{html_escape(row["cam_ccd"])}</td>
<td>{html_escape(row["orbit"])}</td>
<td>{html_escape(row["rms_date"])}</td>
<td>{html_escape(row["ra"])}</td>
<td>{html_escape(row["dec"])}</td>
<td>{html_escape(row.get("l", ""))}</td>
<td>{html_escape(row.get("b", ""))}</td>
<td>{html_escape(row["suspected source type"])}</td>
<td>{html_escape(row["notes"])}</td>
</tr>
<tr><th colspan="2">catalog matches</th><th colspan="3">colors</th><th>Gaia source</th><th>Gaia dist</th><th>PS1 dist</th><th>DES dist</th><th>DECaPS dist</th></tr>
<tr>
<td colspan="2">{html_escape(html_catalog_summary(row))}</td>
<td colspan="3">{html_escape(html_color_summary(row))}</td>
<td>{html_escape(row.get("gaia_source_id", ""))}</td>
<td>{html_escape(row.get("gaia_dist_arcsec", ""))}</td>
<td>{html_escape(row.get("panstarrs_dist_arcsec", ""))}</td>
<td>{html_escape(row.get("des_dist_arcsec", ""))}</td>
<td>{html_escape(row.get("decaps_dist_arcsec", ""))}</td>
</tr>
</table>
{html_escape(row["plot_file"])}
</figcaption>\n"""
            html += f"""  <img src="{html_escape(row["plot_file"])}">\n"""
            html += "</figure>\n"
        html += "</section>\n"

    html += "</body>\n</html>\n"
    with open(outfile, "w", newline="") as f:
        f.write(html)


if __name__ == "__main__":
    main()
