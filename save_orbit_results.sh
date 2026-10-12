#!/usr/bin/env bash
# Save one orbit's routed results before the next orbit's analysis overwrites
# them. analysis.sbatch rebuilds sector<N>/{single,multiple,flag} for the orbit
# it is given, and rms_source.py rewrites each cam/ccd's phot.data, so a
# multi-orbit run would otherwise keep only the last orbit.
#
# Run from the directory holding sector<N>/ (WORK_DIR):
#   bash save_orbit_results.sh <sector> <orbit>
# Copies to sector<N>/results/<orbit>/: single/ multiple/ flag/, each cam/ccd's
# phot.data (as phot.data.camC_ccdD), and the multiple_candidates CSV/HTML (with
# the HTML's image links pointing at the copied plots).

set -eo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: bash $0 <sector> <orbit>"
    exit 1
fi
sector=$1
orbit=$2
sdir="sector${sector}"
out="$sdir/results/$orbit"

rm -rf "$out"
mkdir -p "$out"
for d in single multiple flag; do
    [[ -d "$sdir/$d" ]] && cp -r "$sdir/$d" "$out/"
done
for pf in "$sdir"/cam*_ccd*/phot.data; do
    [[ -s "$pf" ]] && cp "$pf" "$out/phot.data.$(basename "$(dirname "$pf")")"
done
for f in "${MULTIPLE_CSV:-multiple_candidates.csv}" "${MULTIPLE_ENRICHED_CSV:-multiple_candidates_enriched.csv}"; do
    [[ -f "$f" ]] && cp "$f" "$out/"
done
html=${MULTIPLE_HTML:-multiple_candidates.html}
if [[ -f "$html" ]]; then
    sed -E "s#${sdir}/(single|multiple|flag)/#\1/#g" "$html" > "$out/$(basename "$html")"
fi

n=$(find "$out" -name '*.jpg' | wc -l)
echo "Saved $orbit results ($n plots) to $out"
