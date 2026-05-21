#!/usr/bin/env bash
# cleanup_phot.sh
# Interactively removes stale phot.data (one per orbit dir, fast) and
# aggregated lc_* files from the local sector directory.
# The slow per-slice lc_* removal is handled by cleanup_lc.sbatch.
#
# Usage:
#   bash cleanup_phot.sh <sector> <cam> <ccd> <lcdir>

set -eo pipefail

if [[ $# -ne 4 ]]; then
    echo "Usage: $0 <sector> <cam> <ccd> <lcdir>"
    exit 1
fi

sector=$1
cam=$2
ccd=$3
lcdir=$4

sectoruse=$(printf "s%04d" "$sector")
dtarget="${DATA_DIR}/${sectoruse}/cam${cam}-ccd${ccd}"
dhome=$(pwd)
outbase="$dhome/sector${sector}/cam${cam}_ccd${ccd}"

echo "======================================================"
echo " Cleanup for sector=$sector cam=$cam ccd=$ccd lcdir=$lcdir"
echo " DATA_DIR target : $dtarget"
echo " Local output    : $outbase"
echo "======================================================"

if [[ ! -d "$dtarget" ]]; then
    echo "ERROR: Target directory not found: $dtarget"
    exit 1
fi

# --- 1. Remove phot.data from all orbit directories (fast) ---
echo ""
echo "--- Removing phot.data from orbit directories ---"
phot_count=0
for o_dir in "$dtarget"/o??; do
    [[ -d "$o_dir" ]] || continue
    pd="$o_dir/phot.data"
    if [[ -e "$pd" ]]; then
        echo "  Removing: $pd"
        rm -f "$pd"
        (( phot_count++ )) || true
    fi
done
echo "  Removed $phot_count phot.data file(s)."

# --- 2. Remove aggregated lcdir from local sector directory (fast) ---
echo ""
echo "--- Removing aggregated lcdir from local sector directory ---"
for target_dir in "$outbase/$lcdir" "$outbase/bkg_phot/$lcdir"; do
    if [[ -d "$target_dir" ]]; then
        echo "  Removing: $target_dir"
        rm -rf "$target_dir"
    fi
done
echo "  Done."

echo ""
echo "Done. Submit cleanup_lc.sbatch separately for per-slice lc_* removal."
