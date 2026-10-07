#!/usr/bin/env python

import argparse
import csv
import os
import time
import warnings

import numpy as np

try:
    from tqdm import tqdm
except ImportError:
    def tqdm(iterable, **kwargs):
        return iterable


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

RA_FIELDS = ["ra", "RA", "RA_deg", "ra_deg"]
DEC_FIELDS = ["dec", "DEC", "Dec", "DEC_deg", "dec_deg"]

MAG_COLUMNS = {
    "g": ["gMeanPSFMag", "gMeanKronMag", "gPSFMag", "gmag", "g_mag", "g"],
    "r": ["rMeanPSFMag", "rMeanKronMag", "rPSFMag", "rmag", "r_mag", "r"],
    "i": ["iMeanPSFMag", "iMeanKronMag", "iPSFMag", "imag", "i_mag", "i"],
    "z": ["zMeanPSFMag", "zMeanKronMag", "zPSFMag", "zmag", "z_mag", "z"],
    "y": ["yMeanPSFMag", "yMeanKronMag", "yPSFMag", "Ymag", "ymag", "y_mag", "y"],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Query Gaia, Pan-STARRS, DES DR2, and DECaPS around exported "
            "multiple-candidate positions and append match/color columns."
        )
    )
    parser.add_argument(
        "infile",
        nargs="?",
        default="multiple_candidates.csv",
        help="Input CSV from export_multiple_candidates.py. Default: multiple_candidates.csv",
    )
    parser.add_argument(
        "--output",
        default="multiple_candidates_enriched.csv",
        help="Output enriched CSV. Default: multiple_candidates_enriched.csv",
    )
    parser.add_argument(
        "--sleep",
        type=float,
        default=0.4,
        help="Seconds to sleep between candidates. Default: 0.4",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Only process the first N input rows. Useful for testing.",
    )
    parser.add_argument(
        "--radius-arcsec",
        type=float,
        default=4.0,
        help="Gaia/DES search radius in arcsec. Default: 4.",
    )
    parser.add_argument(
        "--panstarrs-radius-arcsec",
        type=float,
        default=2.0,
        help="Pan-STARRS search radius in arcsec. Default: 2.",
    )
    parser.add_argument(
        "--decaps-box-arcmin",
        type=float,
        default=3.0,
        help="Half-width of the DECaPS RA/Dec box search in arcmin. Default: 3.",
    )
    parser.add_argument(
        "--decaps-match-radius-arcsec",
        type=float,
        default=4.0,
        help=(
            "Maximum nearest-DECaPS distance to report as a counterpart match. "
            "The broader DECaPS box count is still recorded. Default: 4."
        ),
    )
    parser.add_argument("--no-gaia", action="store_true", help="Skip Gaia DR3.")
    parser.add_argument("--no-panstarrs", action="store_true", help="Skip Pan-STARRS DR2.")
    parser.add_argument("--no-des", action="store_true", help="Skip DES DR2.")
    parser.add_argument("--no-decaps", action="store_true", help="Skip DECaPS DR2.")
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-query rows even if the output CSV already has catalog data.",
    )
    parser.add_argument(
        "--checkpoint-every",
        type=int,
        default=25,
        help="Write partial output every N processed rows. Default: 25.",
    )
    return parser.parse_args()


def clean_value(value):
    if value is None or np.ma.is_masked(value):
        return None
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        value = value.decode("utf-8", errors="replace")
    if isinstance(value, str):
        text = value.strip()
        if text in ("", "--", "nan", "NaN"):
            return None
        return text
    try:
        if not np.isfinite(value):
            return None
    except TypeError:
        pass
    return value


def as_float(value):
    value = clean_value(value)
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def fmt(value, digits=4):
    value = clean_value(value)
    if value is None:
        return ""
    number = as_float(value)
    if number is not None:
        return f"{number:.{digits}f}"
    return str(value)


def fmt_identifier(value):
    value = clean_value(value)
    if value is None:
        return ""
    try:
        return str(int(value))
    except (TypeError, ValueError, OverflowError):
        return str(value)


def get_row_float(row, names):
    for name in names:
        value = as_float(row.get(name, ""))
        if value is not None:
            return value
    return None


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


def read_csv(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    return fieldnames, rows


def write_csv(path, fieldnames, rows):
    tmp_path = f"{path}.tmp"
    with open(tmp_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp_path, path)


def merge_existing_output(rows, output_path):
    if not os.path.exists(output_path) or os.path.getsize(output_path) == 0:
        return
    _, existing_rows = read_csv(output_path)
    existing_by_key = {row_key(row): row for row in existing_rows}
    for row in rows:
        previous = existing_by_key.get(row_key(row))
        if previous is None:
            continue
        for field in CATALOG_FIELDNAMES:
            if previous.get(field, "") != "":
                row[field] = previous[field]


def has_catalog_data(row):
    return any(row.get(field, "") for field in CATALOG_FIELDNAMES)


def ensure_fields(fieldnames):
    out = list(fieldnames)
    for field in CATALOG_FIELDNAMES:
        if field not in out:
            out.append(field)
    return out


def table_column(table_or_row, names):
    colnames = getattr(table_or_row, "colnames", None)
    if colnames is None and hasattr(table_or_row, "table"):
        colnames = table_or_row.table.colnames
    if colnames is None and hasattr(table_or_row, "_table"):
        colnames = table_or_row._table.colnames
    if colnames is None:
        return None

    exact = set(colnames)
    lower = {name.lower(): name for name in colnames}
    for name in names:
        if name in exact:
            return name
        if name.lower() in lower:
            return lower[name.lower()]
    return None


def record_value(record, names):
    col = table_column(record, names)
    if col is None:
        return None
    try:
        return clean_value(record[col])
    except Exception:
        return None


def color_value(a, b):
    a = as_float(a)
    b = as_float(b)
    if a is None or b is None:
        return ""
    return fmt(a - b)


def set_band_values(outrow, prefix, record):
    values = {}
    for band, names in MAG_COLUMNS.items():
        values[band] = record_value(record, names)
        outrow[f"{prefix}_{band}_mag"] = fmt(values[band])
    outrow[f"{prefix}_g_r"] = color_value(values.get("g"), values.get("r"))
    outrow[f"{prefix}_r_i"] = color_value(values.get("r"), values.get("i"))
    return values


def best_table_index_by_distance(table, sc, dist_names, dist_multiplier, ra_names, dec_names):
    if len(table) == 0:
        return None, None

    dist_col = table_column(table, dist_names)
    if dist_col is not None:
        distances = np.array([as_float(value) for value in table[dist_col]], dtype=float)
        distances *= dist_multiplier
    else:
        ra_col = table_column(table, ra_names)
        dec_col = table_column(table, dec_names)
        if ra_col is None or dec_col is None:
            return 0, None
        from astropy.coordinates import SkyCoord
        import astropy.units as astrou

        coords = SkyCoord(np.array(table[ra_col], dtype=float), np.array(table[dec_col], dtype=float), unit="deg")
        distances = sc.separation(coords).arcsec

    finite = np.isfinite(distances)
    if not finite.any():
        return 0, None
    finite_indices = np.where(finite)[0]
    best_index = int(finite_indices[np.argmin(distances[finite])])
    return best_index, float(distances[best_index])


def add_match(matches, label, distance):
    if distance == "":
        matches.append(label)
    else:
        matches.append(f"{label} {distance} arcsec")


def add_color(colors, label, value):
    if value:
        colors.append(f"{label}={value}")


def query_gaia(outrow, sc, radius_arcsec, matches, colors):
    from astroquery.gaia import Gaia
    import astropy.units as astrou

    Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    job = Gaia.cone_search_async(sc, radius=radius_arcsec * astrou.arcsec)
    table = job.get_results()
    if len(table) == 0:
        return

    best_idx, dist = best_table_index_by_distance(
        table,
        sc,
        ["dist", "DIST"],
        3600.0,
        ["ra", "RA"],
        ["dec", "DEC"],
    )
    if best_idx is None:
        return
    record = table[best_idx]

    outrow["gaia_source_id"] = fmt_identifier(record_value(record, ["source_id"]))
    outrow["gaia_dist_arcsec"] = fmt(dist)
    outrow["gaia_parallax"] = fmt(record_value(record, ["parallax"]))
    g_mag = record_value(record, ["phot_g_mean_mag"])
    bp_mag = record_value(record, ["phot_bp_mean_mag"])
    rp_mag = record_value(record, ["phot_rp_mean_mag"])
    outrow["gaia_g_mag"] = fmt(g_mag)
    outrow["gaia_bp_mag"] = fmt(bp_mag)
    outrow["gaia_rp_mag"] = fmt(rp_mag)
    outrow["gaia_bp_rp"] = color_value(bp_mag, rp_mag)

    add_match(matches, "Gaia", outrow["gaia_dist_arcsec"])
    add_color(colors, "Gaia BP-RP", outrow["gaia_bp_rp"])


def query_panstarrs(outrow, sc, radius_arcsec, matches, colors):
    from astroquery.mast import Catalogs
    import astropy.units as astrou

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Query returned no results")
        table = Catalogs.query_region(
            sc,
            radius=radius_arcsec * astrou.arcsec,
            catalog="PANSTARRS",
            table="stack",
            data_release="dr2",
        )
    if len(table) == 0:
        return

    best_idx, dist = best_table_index_by_distance(
        table,
        sc,
        ["distance"],
        3600.0,
        ["raMean", "raStack", "ra", "RAJ2000"],
        ["decMean", "decStack", "dec", "DEJ2000"],
    )
    if best_idx is None:
        return
    record = table[best_idx]

    outrow["panstarrs_dist_arcsec"] = fmt(dist)
    set_band_values(outrow, "panstarrs", record)

    add_match(matches, "PS1", outrow["panstarrs_dist_arcsec"])
    add_color(colors, "PS1 g-r", outrow["panstarrs_g_r"])
    add_color(colors, "PS1 r-i", outrow["panstarrs_r_i"])


def query_des(outrow, sc, radius_arcsec, matches, colors):
    from astroquery.vizier import Vizier
    import astropy.units as astrou

    vizier = Vizier(columns=["**"], row_limit=20)
    result = vizier.query_region(
        sc,
        radius=radius_arcsec * astrou.arcsec,
        catalog="II/371/des_dr2",
    )
    table = None
    for candidate in result:
        if len(candidate) > 0:
            table = candidate
            break
    if table is None:
        return

    best_idx, dist = best_table_index_by_distance(
        table,
        sc,
        ["_r"],
        1.0,
        ["RAJ2000", "RA_ICRS", "ra"],
        ["DEJ2000", "DE_ICRS", "dec"],
    )
    if best_idx is None:
        return
    record = table[best_idx]

    outrow["des_dist_arcsec"] = fmt(dist)
    set_band_values(outrow, "des", record)

    add_match(matches, "DES", outrow["des_dist_arcsec"])
    add_color(colors, "DES g-r", outrow["des_g_r"])
    add_color(colors, "DES r-i", outrow["des_r_i"])


def in_decaps_region(row):
    gal_l = get_row_float(row, ["l", "gal_l", "galactic_l"])
    gal_b = get_row_float(row, ["b", "gal_b", "galactic_b"])
    if gal_l is None or gal_b is None:
        return True
    # The original snippet used (l > 124.1) and (l < 6.1), which cannot be
    # satisfied. This treats that as a wrap-around Galactic longitude range.
    return abs(gal_b) < 10.1 and (gal_l > 124.1 or gal_l < 6.1)


def query_decaps(outrow, row, sc, service, box_arcmin, match_radius_arcsec, matches, colors):
    if not in_decaps_region(row):
        return

    ra = sc.ra.deg
    dec = sc.dec.deg
    half_width = box_arcmin / 60.0
    query = (
        "SELECT TOP 50 * FROM decaps_dr2.object "
        f"WHERE ra BETWEEN {ra - half_width} AND {ra + half_width} "
        f"AND dec BETWEEN {dec - half_width} AND {dec + half_width}"
    )
    result = service.search(query)
    table = result.to_table()
    outrow["decaps_match_count"] = str(len(table))
    if len(table) == 0:
        return

    best_idx, dist = best_table_index_by_distance(
        table,
        sc,
        ["dist", "distance"],
        3600.0,
        ["ra", "RA"],
        ["dec", "DEC"],
    )
    if best_idx is None:
        return
    record = table[best_idx]

    if dist is None or dist > match_radius_arcsec:
        return

    outrow["decaps_dist_arcsec"] = fmt(dist)
    set_band_values(outrow, "decaps", record)

    add_match(matches, "DECaPS", outrow["decaps_dist_arcsec"])
    add_color(colors, "DECaPS g-r", outrow["decaps_g_r"])
    add_color(colors, "DECaPS r-i", outrow["decaps_r_i"])


def enrich_row(row, args, decaps_service):
    from astropy.coordinates import SkyCoord

    ra = get_row_float(row, RA_FIELDS)
    dec = get_row_float(row, DEC_FIELDS)
    for field in CATALOG_FIELDNAMES:
        row[field] = ""
    if ra is None or dec is None:
        row["catalog_status"] = "missing ra/dec"
        return

    sc = SkyCoord(ra, dec, unit="deg")
    matches = []
    colors = []
    errors = []

    if not args.no_gaia:
        try:
            query_gaia(row, sc, args.radius_arcsec, matches, colors)
        except Exception as exc:
            errors.append(f"Gaia: {exc}")

    if not args.no_panstarrs:
        try:
            query_panstarrs(row, sc, args.panstarrs_radius_arcsec, matches, colors)
        except Exception as exc:
            errors.append(f"Pan-STARRS: {exc}")

    if not args.no_des:
        try:
            query_des(row, sc, args.radius_arcsec, matches, colors)
        except Exception as exc:
            errors.append(f"DES: {exc}")

    if not args.no_decaps and decaps_service is not None:
        try:
            query_decaps(
                row,
                row,
                sc,
                decaps_service,
                args.decaps_box_arcmin,
                args.decaps_match_radius_arcsec,
                matches,
                colors,
            )
        except Exception as exc:
            errors.append(f"DECaPS: {exc}")

    row["catalog_matches"] = "; ".join(matches)
    row["color_summary"] = "; ".join(colors)
    row["catalog_status"] = "ok" if not errors else "WARNING: " + " | ".join(errors)


def main():
    args = parse_args()
    fieldnames, rows = read_csv(args.infile)
    fieldnames = ensure_fields(fieldnames)
    merge_existing_output(rows, args.output)

    decaps_service = None
    if not args.no_decaps:
        try:
            import pyvo as vo

            decaps_service = vo.dal.TAPService("https://datalab.noirlab.edu/tap")
        except Exception as exc:
            print(f"WARNING: DECaPS disabled because pyvo/TAP setup failed: {exc}")
            args.no_decaps = True

    if args.limit is None:
        process_indices = list(range(len(rows)))
    else:
        process_indices = list(range(min(args.limit, len(rows))))

    processed = 0
    skipped = 0
    for idx in tqdm(process_indices, desc="catalog matches"):
        row = rows[idx]
        if has_catalog_data(row) and not args.force:
            skipped += 1
            continue

        enrich_row(row, args, decaps_service)
        processed += 1
        if args.sleep > 0:
            time.sleep(args.sleep)
        if args.checkpoint_every > 0 and processed % args.checkpoint_every == 0:
            write_csv(args.output, fieldnames, rows)

    write_csv(args.output, fieldnames, rows)
    rows_with_matches = sum(1 for row in rows if row.get("catalog_matches", ""))
    rows_with_gaia = sum(1 for row in rows if row.get("gaia_source_id", ""))
    rows_with_panstarrs = sum(1 for row in rows if row.get("panstarrs_dist_arcsec", ""))
    rows_with_des = sum(1 for row in rows if row.get("des_dist_arcsec", ""))
    rows_with_decaps = sum(1 for row in rows if row.get("decaps_dist_arcsec", ""))
    rows_missing_radec = sum(1 for row in rows if row.get("catalog_status", "") == "missing ra/dec")
    print(f"Wrote {len(rows)} rows to {args.output}")
    print(f"Processed {processed} rows; skipped {skipped} rows with existing catalog data")
    print(
        "Catalog match summary: "
        f"any={rows_with_matches}, Gaia={rows_with_gaia}, "
        f"Pan-STARRS={rows_with_panstarrs}, DES={rows_with_des}, "
        f"DECaPS={rows_with_decaps}, missing_ra_dec={rows_missing_radec}"
    )


if __name__ == "__main__":
    main()
