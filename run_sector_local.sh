#!/usr/bin/env bash
# Whole-sector local run with bounded disk: for each half-orbit in turn,
#   1) upstream/tica_upstream.sh   download -> subtract -> tarballs (raw FFIs deleted)
#   2) orbit_pipeline.sbatch       RMS -> candidates -> light curves -> plots
#   3) purge                       that orbit's tarballs are deleted (PURGE_IMAGES=1)
# so the disk only ever holds one half-orbit of difference images (~18 GB per
# CCD at 10-min cadence, ~50 GB per CCD at 200-s cadence).
#
# Usage:
#   bash run_sector_local.sh <sector> [cam|all] [ccd|all]
# Examples:
#   bash run_sector_local.sh 107 3 3        # one CCD, all four half-orbits
#   bash run_sector_local.sh 107            # every CCD (needs ~16x the disk per orbit)
#
# Rerunning resumes: finished slices and purged orbits are skipped. It stops if
# an orbit's upstream step is incomplete, so a partial orbit is never purged.

set -eo pipefail

if [[ $# -lt 1 || $# -gt 3 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: bash $0 <sector> [cam|all] [ccd|all]"
    exit 1
fi
sector=$1
cam_arg=${2:-all}
ccd_arg=${3:-all}

code_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$code_dir/setup_vars_local.sh"

cams=$([[ "$cam_arg" == "all" ]] && echo 1 2 3 4 || echo "$cam_arg")
ccds=$([[ "$ccd_arg" == "all" ]] && echo 1 2 3 4 || echo "$ccd_arg")
sectoruse=$(printf "s%04d" "$sector")

# Half-orbits come from the manifests (listing only; nothing is downloaded yet).
orbits=()
for cam in $cams; do
    for ccd in $ccds; do
        python "$code_dir/upstream/tica_manifest.py" --sector "$sector" --cam "$cam" --ccd "$ccd" > /dev/null \
            || { echo "WARNING: no TICA listing for cam${cam}-ccd${ccd}"; continue; }
        mapfile -t -O "${#orbits[@]}" orbits < <(awk 'NR>1 {print $4}' \
            "$DATA_DIR/$sectoruse/cam${cam}-ccd${ccd}/manifest.tsv" | sort -u)
    done
done
mapfile -t orbits < <(printf '%s\n' "${orbits[@]}" | sort -u)
if [[ ${#orbits[@]} -eq 0 ]]; then
    echo "ERROR: no TICA frames found for $sectoruse cams=[$cams] ccds=[$ccds]"
    exit 1
fi

echo "======================================================"
echo " Local sector run: $sectoruse  cams: $cams  ccds: $ccds"
echo " Half-orbits: ${orbits[*]}   PURGE_IMAGES=$PURGE_IMAGES"
echo "======================================================"

for orbit in "${orbits[@]}"; do
    echo ""
    echo "################ $sectoruse $orbit: upstream ################"
    if ! bash "$code_dir/upstream/tica_upstream.sh" "$sector" "$cam_arg" "$ccd_arg" "$orbit"; then
        echo "STOPPED: upstream for $orbit is incomplete (see the messages above)."
        echo "Rerun the same command to resume; finished slices are skipped."
        exit 1
    fi
    echo ""
    echo "################ $sectoruse $orbit: orbit pipeline ################"
    (cd "$WORK_DIR" && bash "$code_dir/orbit_pipeline.sbatch" "$sector" "$orbit")
done

echo ""
echo "Sector $sectoruse done. Plots: $WORK_DIR/sector${sector}/{single,multiple,flag}/"
