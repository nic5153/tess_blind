#!/usr/bin/env bash
# Local upstream processing: TICA download -> registration -> image
# subtraction -> HPCC-style products under $DATA_DIR (ref.fits, ref_epsf.fits,
# oXX/rms_list_tica, oXX/sliceNNNN/{dates,images.tar}). Run this before
# orbit_pipeline for the same sector/orbit.
#
# Usage:
#   bash upstream/tica_upstream.sh <sector> [cam|all] [ccd|all] [orbit|all]
#
# cam, ccd and orbit default to all, so "bash upstream/tica_upstream.sh 101"
# processes every cam/ccd/half-orbit of the sector.
#
# Examples:
#   bash upstream/tica_upstream.sh 40 3 1 o1a          # one CCD, one half-orbit
#   bash upstream/tica_upstream.sh 40 all all o1a      # whole focal plane, o1a
#   bash upstream/tica_upstream.sh 40                  # everything in sector 40
#   COMPRESS_Q=4 bash upstream/tica_upstream.sh 40 3 1 o1a   # fpack q=4 test tree
#   MAX_SLICES=3 bash upstream/tica_upstream.sh 40 3 1 o1a   # quick test
#   FORCE=1 ...                                             # redo finished slices
#
# Reruns skip finished slices (.done markers). The reference is built once per
# sector/CCD from REF_ORBIT (default: the first half-orbit, as in the MIT
# TESS-transients convention).
#
# Disk guard: before each half-orbit, the frames still to process (~18 MB of
# conv_ image each) plus MIN_FREE_GB (default 20) must fit on the DATA_DIR
# drive, otherwise the run stops. A whole sector is ~70 GB per CCD at 10-min
# cadence and ~150 GB per CCD at 200-s cadence, so for all 16 CCDs run
# orbit_pipeline and clear finished orbits between upstream runs.

set -eo pipefail

if [[ $# -lt 1 || $# -gt 4 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "Usage: bash $0 <sector> [cam|all] [ccd|all] [orbit|all]   (defaults: all)"
    exit 1
fi

sector=$1
cam_arg=${2:-all}
ccd_arg=${3:-all}
orbit_arg=${4:-all}

if [[ "$orbit_arg" != "all" && ! "$orbit_arg" =~ ^o[0-9][ab]$ ]]; then
    echo "ERROR: orbit must look like o1a, o1b, o2a, o2b (or all)"
    exit 1
fi

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# setup_vars_hpcc.sh exports SETUP_VARS (and TESS_CODE_DIR); default: local.
source "${SETUP_VARS:-$script_dir/../setup_vars_local.sh}"
[[ -n "${TESS_CODE_DIR:-}" && -d "$TESS_CODE_DIR/upstream" ]] && script_dir="$TESS_CODE_DIR/upstream"

cams=$([[ "$cam_arg" == "all" ]] && echo 1 2 3 4 || echo "$cam_arg")
ccds=$([[ "$ccd_arg" == "all" ]] && echo 1 2 3 4 || echo "$ccd_arg")
sectoruse=$(printf "s%04d" "$sector")
timestamp=$(date +%Y%m%d_%H%M%S)
failures=0

echo "======================================================"
echo " Local TICA upstream processing"
echo " Sector : $sectoruse   cams: $cams   ccds: $ccds   orbit: $orbit_arg"
echo " DATA_DIR : $DATA_DIR"
echo " RAW_CACHE: $RAW_CACHE"
echo " Workers  : $NPROC   SUB_MODE=${SUB_MODE:-reference}   COMPRESS_Q=${COMPRESS_Q:-off}"
echo "======================================================"

for cam in $cams; do
    for ccd in $ccds; do
        tag="${sectoruse}_cam${cam}-ccd${ccd}"
        log="$LOG_DIR/upstream_${tag}_${orbit_arg}_${timestamp}.log"
        echo ""
        echo "--- cam${cam}-ccd${ccd}  (log: $log)"

        if ! python "$script_dir/tica_manifest.py" --sector "$sector" --cam "$cam" --ccd "$ccd" 2>&1 | tee -a "$log"; then
            echo "WARNING: no TICA listing for cam${cam}-ccd${ccd}; skipping" | tee -a "$log"
            failures=$((failures + 1))
            continue
        fi
        manifest="$DATA_DIR/$sectoruse/cam${cam}-ccd${ccd}/manifest.tsv"

        mapfile -t orbits < <(awk 'NR>1 {print $4}' "$manifest" | sort -u)
        if [[ "$orbit_arg" != "all" ]]; then
            if [[ ! " ${orbits[*]} " =~ " $orbit_arg " ]]; then
                echo "WARNING: $orbit_arg not in manifest (have: ${orbits[*]}); skipping" | tee -a "$log"
                continue
            fi
            orbits=("$orbit_arg")
        fi

        ref_orbit=${REF_ORBIT:-$(awk 'NR>1 {print $4}' "$manifest" | sort -u | head -1)}
        if ! python "$script_dir/tica_make_ref.py" --sector "$sector" --cam "$cam" --ccd "$ccd" \
                --orbit "$ref_orbit" 2>&1 | tee -a "$log" \
            || ! python "$script_dir/tica_make_epsf.py" --sector "$sector" --cam "$cam" --ccd "$ccd" \
                2>&1 | tee -a "$log"; then
            echo "WARNING: reference/ePSF failed for cam${cam}-ccd${ccd}; skipping" | tee -a "$log"
            failures=$((failures + 1))
            continue
        fi

        for orbit in "${orbits[@]}"; do
            if [[ -e "$DATA_DIR/$sectoruse/cam${cam}-ccd${ccd}/$orbit/.purged" ]]; then
                echo "  $orbit: already processed by orbit_pipeline and purged; skipping" | tee -a "$log"
                continue
            fi
            mapfile -t slices < <(awk -v o="$orbit" 'NR>1 && $4==o {print $5}' "$manifest" | sort -u)
            if [[ -n "${MAX_SLICES:-}" ]]; then
                slices=("${slices[@]:0:$MAX_SLICES}")
            fi

            # Disk guard: frames in unfinished slices x ~18 MB, plus a margin.
            orbit_path="$DATA_DIR/$sectoruse/cam${cam}-ccd${ccd}/$orbit"
            todo_frames=0
            while read -r sl n; do
                if [[ " ${slices[*]} " == *" $sl "* && ! -e "$orbit_path/$sl/.done" ]]; then
                    todo_frames=$((todo_frames + n))
                fi
            done < <(awk -v o="$orbit" 'NR>1 && $4==o {n[$5]++} END {for (s in n) print s, n[s]}' "$manifest")
            need_gb=$(( todo_frames * 18 / 1024 + ${MIN_FREE_GB:-20} ))
            free_gb=$(( $(df -Pk "$DATA_DIR" | awk 'NR==2 {print $4}') / 1024 / 1024 ))
            if (( todo_frames > 0 && free_gb < need_gb )); then
                echo "ERROR: $orbit needs ~${need_gb} GB ($todo_frames frames + ${MIN_FREE_GB:-20} GB margin)" \
                     "but only ${free_gb} GB is free on $DATA_DIR. Run orbit_pipeline on finished" \
                     "orbits and clear their slice tarballs, then rerun." | tee -a "$log"
                exit 1
            fi
            echo "  $orbit: ${#slices[@]} slices ($todo_frames frames to do, ~$((todo_frames * 18 / 1024)) GB;" \
                 "${free_gb} GB free) with $NPROC workers" | tee -a "$log"

            force_flag=$([[ "${FORCE:-0}" == "1" ]] && echo "--force" || true)
            printf '%s\n' "${slices[@]}" | xargs -P "$NPROC" -I{} bash -c '
                python "$0/tica_process_slice.py" --sector "$1" --cam "$2" --ccd "$3" \
                    --orbit "$4" --slice "{}" $5 || echo "FAILED: $4/{}"
            ' "$script_dir" "$sector" "$cam" "$ccd" "$orbit" "$force_flag" 2>&1 | tee -a "$log" | grep -E "Done|FAILED|SKIP|ERROR|Traceback" || true

            if ! python "$script_dir/tica_finalize_orbit.py" --sector "$sector" --cam "$cam" \
                    --ccd "$ccd" --orbit "$orbit" 2>&1 | tee -a "$log"; then
                failures=$((failures + 1))
            fi
        done
    done
done

echo ""
if [[ $failures -gt 0 ]]; then
    echo "Finished with $failures incomplete orbit(s); rerun the same command to retry missing slices."
    exit 1
fi
echo "Upstream products ready under $DATA_DIR/$sectoruse"
