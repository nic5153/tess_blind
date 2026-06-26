#!/usr/bin/env bash
#SBATCH --job-name="source"
#SBATCH --output=/lustre/work/nimcclur/TESS/photometry/logs/%x.o%A_%a
#SBATCH --error=/lustre/work/nimcclur/TESS/photometry/logs/%x.e%A_%a
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=0-4:00:00
#SBATCH --array=0-15%8

set -eo pipefail

if [[ $# -ne 1 && $# -ne 3 && $# -ne 4 ]]; then
    echo "Usage:"
    echo "  sbatch $0 <sector>                       # whole sector; use array task IDs"
    echo "  sbatch $0 <sector> <cam> <ccd>           # one cam/ccd; submit with --array=0-0"
    echo "  sbatch $0 <sector> <cam> <ccd> <orbit>   # one cam/ccd/orbit; submit with --array=0-0"
    exit 1
fi

sector=$1
requested_cam=${2:-}
requested_ccd=${3:-}
requested_orbit=${4:-}
sector_dir="sector${sector}"

if [[ -d "$sector_dir" ]]; then
    cd "$sector_dir"
elif [[ "$(basename "$PWD")" == "$sector_dir" ]]; then
    :
else
    echo "ERROR: run from /home/nimcclur/work/TESS/photometry or sector${sector}"
    exit 1
fi

mkdir -p /lustre/work/nimcclur/TESS/photometry/logs

if [[ ! -f do_photutils_extract.py ]]; then
    echo "ERROR: missing do_photutils_extract.py in $(pwd)"
    exit 1
fi

if [[ ! -f rms_source.py ]]; then
    echo "ERROR: missing rms_source.py in $(pwd)"
    exit 1
fi

if [[ -n "$requested_orbit" && -z "$requested_cam" ]]; then
    echo "ERROR: orbit mode requires cam and ccd"
    exit 1
fi

if [[ -n "$requested_orbit" && ! "$requested_orbit" =~ ^o[0-9][ab]$ ]]; then
    echo "ERROR: orbit must look like o1a, o1b, o2a, o2b"
    exit 1
fi

if [[ -n "$requested_cam" && -n "$requested_ccd" ]]; then
    camccd_dir="cam${requested_cam}_ccd${requested_ccd}"

    if [[ ! -d "$camccd_dir" ]]; then
        echo "ERROR: requested cam/ccd directory not found: $camccd_dir"
        exit 1
    fi

    if [[ "${SLURM_ARRAY_TASK_ID:-0}" -ne 0 ]]; then
        echo "Single target mode should be submitted with --array=0-0; exiting duplicate task ${SLURM_ARRAY_TASK_ID:-0}."
        exit 0
    fi
else
    mapfile -t camccd_dirs < <(find . -maxdepth 1 -type d -name 'cam*_ccd*' -printf '%f\n' | sort)

    if [[ ${#camccd_dirs[@]} -eq 0 ]]; then
        echo "ERROR: no cam*_ccd* dirs found in $(pwd)"
        exit 1
    fi

    if [[ "$SLURM_ARRAY_TASK_ID" -ge "${#camccd_dirs[@]}" ]]; then
        echo "Task $SLURM_ARRAY_TASK_ID has no cam/ccd directory; exiting."
        exit 0
    fi

    camccd_dir="${camccd_dirs[$SLURM_ARRAY_TASK_ID]}"
fi

if [[ "$camccd_dir" =~ cam([0-9]+)_ccd([0-9]+) ]]; then
    cam="${BASH_REMATCH[1]}"
    ccd="${BASH_REMATCH[2]}"
else
    echo "ERROR: cannot parse cam/ccd from $camccd_dir"
    exit 1
fi

if [[ -n "$requested_orbit" ]]; then
    orbit_dirs=("$camccd_dir/$requested_orbit")
    if [[ ! -d "${orbit_dirs[0]}" ]]; then
        echo "ERROR: requested orbit directory not found: ${orbit_dirs[0]}"
        exit 1
    fi
else
    orbit_dirs=("$camccd_dir"/o??)
fi

echo "======================================================"
echo "Task ${SLURM_ARRAY_TASK_ID:-0} processing $camccd_dir"
echo "sector=$sector cam=$cam ccd=$ccd orbit=${requested_orbit:-all}"

for orbit_dir in "${orbit_dirs[@]}"; do
    [[ -d "$orbit_dir" ]] || continue
    echo "--- Running do_photutils_extract.py in $orbit_dir"

    (
        cd "$orbit_dir"
        shopt -s nullglob
        for rms_fits in rms*.fits; do
            outcat="outcat${rms_fits%.fits}.tsv"
            if [[ -s "$outcat" ]]; then
                echo "SKIP existing $PWD/$outcat"
                continue
            fi

            echo "python ../../do_photutils_extract.py $rms_fits"
            python ../../do_photutils_extract.py "$rms_fits"
        done
    )
done

echo "--- Running rms_source.py for $camccd_dir"
(
    cd "$camccd_dir"
    shopt -s nullglob

    if [[ -n "$requested_orbit" ]]; then
        outcats=("$requested_orbit"/outcatrms*.tsv)
    else
        outcats=(o??/outcatrms*.tsv)
    fi

    if [[ ${#outcats[@]} -eq 0 ]]; then
        echo "ERROR: no outcatrms*.tsv files found for $camccd_dir"
        exit 1
    fi

    echo "python ../rms_source.py --sector $sector --cam $cam --ccd $ccd --rms ${outcats[*]}"
    python ../rms_source.py --sector "$sector" --cam "$cam" --ccd "$ccd" --rms "${outcats[@]}"
)

echo "Done with $camccd_dir"
