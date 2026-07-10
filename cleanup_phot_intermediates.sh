#!/usr/bin/env bash
#
# Remove loose photometry intermediate products left in slice directories.
# This mirrors the cleanup performed inside do_phot_em2.sbatch:
#   rm -f -- {interp,conv,bkg}*fits kt_*fits kc_*fits
#
# Default mode is a dry run. Add --delete to actually remove files.

set -eo pipefail

usage() {
    cat <<'EOF'
Usage:
  bash cleanup_phot_intermediates.sh <sector> [cam] [ccd] [orbit] [options]

Examples:
  # Dry-run every matching file from sector 102 owned by you
  bash cleanup_phot_intermediates.sh 102

  # Dry-run one cam/ccd/orbit
  bash cleanup_phot_intermediates.sh 102 1 4 o2a

  # Actually delete one cam/ccd/orbit
  bash cleanup_phot_intermediates.sh 102 1 4 o2a --delete

  # Use an explicit data root instead of DATA_DIR
  bash cleanup_phot_intermediates.sh 102 --data-dir /lustre/research/mfausnau/data/tica

Options:
  --delete       Actually remove files. Without this, the script only prints.
  --all-owners   Include files not owned by the current user.
  --data-dir DIR Data root. Defaults to DATA_DIR, then /lustre/research/mfausnau/data/tica.
  -h, --help     Show this help.

Targets:
  Files or symlinks under slice directories matching:
    interp*fits, conv*fits, bkg*fits, kt_*fits, kc_*fits
EOF
}

delete=0
owner_only=1
data_dir="${DATA_DIR:-/lustre/research/mfausnau/data/tica}"
positional=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --delete)
            delete=1
            shift
            ;;
        --all-owners)
            owner_only=0
            shift
            ;;
        --data-dir)
            if [[ $# -lt 2 ]]; then
                echo "ERROR: --data-dir requires a path" >&2
                exit 2
            fi
            data_dir=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        -*)
            echo "ERROR: unknown option: $1" >&2
            usage
            exit 2
            ;;
        *)
            positional+=("$1")
            shift
            ;;
    esac
done

if [[ ${#positional[@]} -lt 1 || ${#positional[@]} -gt 4 ]]; then
    usage
    exit 2
fi

sector=${positional[0]}
cam_filter=${positional[1]:-*}
ccd_filter=${positional[2]:-*}
orbit_filter=${positional[3]:-*}

if [[ ! "$sector" =~ ^[0-9]+$ ]]; then
    echo "ERROR: sector must be an integer, got: $sector" >&2
    exit 2
fi
if [[ "$cam_filter" != "*" && ! "$cam_filter" =~ ^[1-4]$ ]]; then
    echo "ERROR: cam must be 1-4, got: $cam_filter" >&2
    exit 2
fi
if [[ "$ccd_filter" != "*" && ! "$ccd_filter" =~ ^[1-4]$ ]]; then
    echo "ERROR: ccd must be 1-4, got: $ccd_filter" >&2
    exit 2
fi
if [[ "$orbit_filter" != "*" && ! "$orbit_filter" =~ ^o[0-9][ab]$ ]]; then
    echo "ERROR: orbit must look like o1a, o1b, o2a, got: $orbit_filter" >&2
    exit 2
fi

sectoruse=$(printf "s%04d" "$sector")
base="${data_dir%/}/${sectoruse}"

if [[ ! -d "$base" ]]; then
    echo "ERROR: sector directory does not exist: $base" >&2
    exit 1
fi

resolved_base=$(readlink -f "$base")
if [[ -z "$resolved_base" || "$resolved_base" == "/" || "$resolved_base" != */tica/"$sectoruse" ]]; then
    echo "ERROR: refusing to operate on unexpected path: $resolved_base" >&2
    echo "Expected a sector path ending in /tica/$sectoruse" >&2
    exit 1
fi

tmpfile=$(mktemp)
trap 'rm -f "$tmpfile"' EXIT

echo "======================================================"
echo " Photometry intermediate cleanup"
echo " Sector    : $sectoruse"
echo " Data root : ${data_dir%/}"
echo " Cam filter: $cam_filter"
echo " CCD filter: $ccd_filter"
echo " Orbit     : $orbit_filter"
echo " Owner     : $([[ "$owner_only" == "1" ]] && echo "$USER only" || echo "all owners")"
echo " Mode      : $([[ "$delete" == "1" ]] && echo DELETE || echo DRY-RUN)"
echo "======================================================"

shopt -s nullglob
roots=()
for cam_dir in "$base"/cam${cam_filter}-ccd${ccd_filter}; do
    [[ -d "$cam_dir" ]] || continue
    if [[ "$orbit_filter" == "*" ]]; then
        for orbit_dir in "$cam_dir"/o[0-9][ab]; do
            [[ -d "$orbit_dir" ]] || continue
            roots+=("$orbit_dir")
        done
    else
        orbit_dir="$cam_dir/$orbit_filter"
        [[ -d "$orbit_dir" ]] || continue
        roots+=("$orbit_dir")
    fi
done

if [[ ${#roots[@]} -eq 0 ]]; then
    echo "No matching cam/ccd/orbit directories found."
    exit 0
fi

for root in "${roots[@]}"; do
    if [[ "$owner_only" == "1" ]]; then
        find "$root" \
            -path '*/slice*/*' \
            \( -type f -o -type l \) \
            -user "$USER" \
            \( -name 'interp*fits' -o -name 'conv*fits' -o -name 'bkg*fits' -o -name 'kt_*fits' -o -name 'kc_*fits' \) \
            -print0 >> "$tmpfile"
    else
        find "$root" \
            -path '*/slice*/*' \
            \( -type f -o -type l \) \
            \( -name 'interp*fits' -o -name 'conv*fits' -o -name 'bkg*fits' -o -name 'kt_*fits' -o -name 'kc_*fits' \) \
            -print0 >> "$tmpfile"
    fi
done

count=0
bytes=0
while IFS= read -r -d '' file; do
    count=$((count + 1))
    size=$(stat -c '%s' "$file" 2>/dev/null || echo 0)
    bytes=$((bytes + size))
done < "$tmpfile"

human_bytes=$(numfmt --to=iec --suffix=B "$bytes" 2>/dev/null || echo "${bytes} bytes")

echo "Matched files: $count"
echo "Matched size : $human_bytes"

if [[ "$count" -eq 0 ]]; then
    echo "Nothing to clean."
    exit 0
fi

if [[ "$delete" != "1" ]]; then
    echo ""
    echo "DRY RUN: files that would be removed:"
    tr '\0' '\n' < "$tmpfile"
    echo ""
    echo "Re-run with --delete to remove these files."
    exit 0
fi

echo ""
echo "Deleting matched files..."
while IFS= read -r -d '' file; do
    rm -f -- "$file"
done < "$tmpfile"

echo "Deleted $count files ($human_bytes)."
