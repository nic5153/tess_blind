#!/usr/bin/env bash
set -euo pipefail

remote_user="${HPCC_USER:-nimcclur}"
remote_host="${HPCC_HOST:-login.hpcc.ttu.edu}"
remote_root="${HPCC_PHOTOMETRY_ROOT:-/work/TESS/photometry}"
outdir="."
csv_name="multiple_candidates.csv"
html_name="multiple_candidates.html"

usage() {
  cat <<EOF
Usage: bash download_multiple_candidates.sh [options]

Downloads the multiple-candidate CSV and HTML generated on HPCC.

Options:
  --user USER          HPCC username. Default: HPCC_USER or nimcclur
  --host HOST          HPCC login host. Default: HPCC_HOST or login.hpcc.ttu.edu
  --remote-root PATH   Remote photometry root. Default: HPCC_PHOTOMETRY_ROOT or /work/TESS/photometry
  --outdir PATH        Local output directory. Default: current directory
  --csv NAME           CSV filename. Default: multiple_candidates.csv
  --html NAME          HTML filename. Default: multiple_candidates.html
  -h, --help           Show this help message
EOF
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --user)
      remote_user="$2"
      shift 2
      ;;
    --host)
      remote_host="$2"
      shift 2
      ;;
    --remote-root)
      remote_root="$2"
      shift 2
      ;;
    --outdir)
      outdir="$2"
      shift 2
      ;;
    --csv)
      csv_name="$2"
      shift 2
      ;;
    --html)
      html_name="$2"
      shift 2
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      echo "Unknown option: $1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

mkdir -p "$outdir"
remote="${remote_user}@${remote_host}"
remote_root="${remote_root%/}"
remote_csv="${remote_root}/${csv_name}"
remote_html="${remote_root}/${html_name}"

echo "Downloading from ${remote}:${remote_root}"

if ! ssh "$remote" "test -f '$remote_csv' && test -f '$remote_html'"; then
  cat >&2 <<EOF
Missing one or both remote exporter outputs:
  ${remote}:${remote_csv}
  ${remote}:${remote_html}

Generate them on HPCC first:
  ssh ${remote}
  cd ${remote_root}
  python export_multiple_candidates.py

If your photometry repo is in a different HPCC path, rerun this script with:
  bash download_multiple_candidates.sh --remote-root /path/to/photometry
EOF
  exit 1
fi

scp "${remote}:${remote_csv}" "${outdir}/${csv_name}"
scp "${remote}:${remote_html}" "${outdir}/${html_name}"

echo
echo "Downloaded:"
echo "  ${outdir}/${csv_name}"
echo "  ${outdir}/${html_name}"
