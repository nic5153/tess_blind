#!/usr/bin/env bash
# Re-make the candidate list, plots and exports for orbits that were already
# processed (their light curves exist), without photometry or image products:
#   source.sbatch (re-extracts missing outcatrms_day_*.tsv from the saved daily
#                  RMS images, then rms_source.py) -> phot.data
#   analysis.sbatch -> plots routed to single/multiple/flag
#   export_multiple_candidates.py -> CSV/HTML
#   save_orbit_results.sh -> sector<N>/results/<orbit>/
#
# Usage (local):  bash reanalyze_orbit.sh <sector> <orbit> [orbit ...]
# Usage (HPCC):   source setup_vars_hpcc.sh   (or ~/setup_vars.sh with SETUP_VARS set)
#                 sbatch -p nocona -c 8 -t 1-00:00:00 --wrap "bash $TESS_CODE_DIR/reanalyze_orbit.sh 101 o1a o1b"
# If the orbit's difference images were purged, plots lack the detection-image
# cutout and Gaussian-fit panels; everything else is unchanged.

set -eo pipefail

if [[ $# -lt 2 ]]; then
    echo "Usage: bash $0 <sector> <orbit> [orbit ...]"
    exit 1
fi
sector=$1
shift

code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SETUP_VARS:-$code_dir/setup_vars_local.sh}"
code_dir=${TESS_CODE_DIR:-$code_dir}
cd "$WORK_DIR"
sdir="sector${sector}"
export ANALYSIS_N=${ANALYSIS_N:-$NPROC}

for orbit in "$@"; do
    echo "=== $sdir $orbit ==="
    for camdir in "$sdir"/cam*_ccd*; do
        compgen -G "$camdir/$orbit/rms_day_*.fits" > /dev/null || continue
        [[ "$(basename "$camdir")" =~ cam([0-9])_ccd([0-9]) ]] || continue
        cam=${BASH_REMATCH[1]}
        ccd=${BASH_REMATCH[2]}
        if bash "$code_dir/source.sbatch" "$sector" "$cam" "$ccd" "$orbit" \
                > "$LOG_DIR/reanalyze_source_${sdir}_cam${cam}_ccd${ccd}_${orbit}.log" 2>&1; then
            echo "  cam${cam}_ccd${ccd}: phot.data rebuilt ($(wc -l < "$camdir/phot.data") candidates)"
        else
            echo "  WARNING: source.sbatch failed for cam${cam}_ccd${ccd} (see $LOG_DIR)"
        fi
    done
    bash "$code_dir/analysis.sbatch" "$sector" "$orbit" > "$LOG_DIR/reanalyze_${sdir}_${orbit}.log" 2>&1 \
        || echo "WARNING: analysis.sbatch reported an error (see $LOG_DIR/reanalyze_${sdir}_${orbit}.log)"
    python "$code_dir/export_multiple_candidates.py" --root "." \
        --outfile "${MULTIPLE_CSV:-multiple_candidates.csv}" \
        --html "${MULTIPLE_HTML:-multiple_candidates.html}" \
        --enriched-csv "${MULTIPLE_ENRICHED_CSV:-multiple_candidates_enriched.csv}" \
        || echo "WARNING: export_multiple_candidates.py failed"
    bash "$code_dir/save_orbit_results.sh" "$sector" "$orbit"
done
