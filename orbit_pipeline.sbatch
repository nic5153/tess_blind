#!/usr/bin/env bash

#SBATCH --job-name=orbit_tess_pipeline
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00

set -eo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
    echo "Usage: sbatch $0 <sector> <orbit> [lcdir]"
    echo "Example: sbatch $0 105 o1a lc"
    exit 1
fi

sector=$1
orbit=$2
lcdir=${3:-lc}

if [[ ! "$orbit" =~ ^o[0-9][ab]$ ]]; then
    echo "ERROR: orbit must look like o1a, o1b, o2a, o2b"
    exit 1
fi

if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/do_photutils_extract.py" ]]; then
    script_dir="$SLURM_SUBMIT_DIR"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source /home/${USER}/setup_vars.sh
logdir="${LOG_DIR:-/lustre/work/nimcclur/TESS/photometry/logs}"
mkdir -p "$logdir"

sectoruse=$(printf "s%04d" "$sector")
dtarget_base="${DATA_DIR}/${sectoruse}"

wait_for_jobs() {
    local label=$1
    shift
    local jobs=("$@")

    if [[ ${#jobs[@]} -eq 0 ]]; then
        echo "No jobs to wait for: $label"
        return 0
    fi

    echo "Waiting for $label: ${jobs[*]}"
    while true; do
        local running=0
        for job in "${jobs[@]}"; do
            if squeue -h -j "$job" | grep -q .; then
                running=1
                break
            fi
        done
        [[ "$running" -eq 0 ]] && break
        sleep 60
    done
    echo "Finished waiting for $label"
}

echo "======================================================"
echo " Incremental TESS blind orbit pipeline"
echo " Sector : $sectoruse"
echo " Orbit  : $orbit"
echo " lcdir  : $lcdir"
echo " CWD    : $(pwd)"
echo " Data   : $dtarget_base"
echo " Scripts: $script_dir"
echo " Logs   : $logdir"
echo "======================================================"

if [[ ! -d "$dtarget_base" ]]; then
    echo "ERROR: sector data directory not found: $dtarget_base"
    exit 1
fi

mkdir -p "sector${sector}"
cp "$script_dir/do_photutils_extract.py" "sector${sector}/do_photutils_extract.py"
cp "$script_dir/rms_source.py" "sector${sector}/rms_source.py"

available_camccds=()
for cam in {1..4}; do
    for ccd in {1..4}; do
        if [[ -d "$dtarget_base/cam${cam}-ccd${ccd}/$orbit" ]]; then
            mkdir -p "sector${sector}/cam${cam}_ccd${ccd}/$orbit"
            available_camccds+=("${cam}:${ccd}")
        else
            echo "Skipping missing cam${cam}_ccd${ccd}/$orbit"
        fi
    done
done

if [[ ${#available_camccds[@]} -eq 0 ]]; then
    echo "ERROR: no cam/ccd directories found for $sectoruse $orbit"
    exit 1
fi

echo ""
echo "[1/5] Generating/updating daily RMS files for available cam/ccd directories..."
rms_jobs=()
for camccd in "${available_camccds[@]}"; do
    cam=${camccd%%:*}
    ccd=${camccd##*:}
    rms_job=$(sbatch --parsable "$script_dir/make_rms_sbatch.sh" "$sector" "$cam" "$ccd" "$orbit")
    echo "  daily_rms cam${cam}_ccd${ccd}: $rms_job"
    rms_jobs+=("$rms_job")
done
wait_for_jobs "daily RMS generation" "${rms_jobs[@]}"

echo ""
echo "[2/5] Refreshing cumulative source lists from all available RMS catalogs..."
source_jobs=()
for camccd in "${available_camccds[@]}"; do
    cam=${camccd%%:*}
    ccd=${camccd##*:}
    source_job=$(sbatch --parsable --array=0-0 "$script_dir/source.sbatch" "$sector" "$cam" "$ccd")
    echo "  source cam${cam}_ccd${ccd}: $source_job"
    source_jobs+=("$source_job")
done
wait_for_jobs "source extraction" "${source_jobs[@]}"

echo ""
echo "[3/5] Running photometry only for $orbit, using cumulative phot.data files..."
copy_jobs=()
do_jobs=()
for camccd in "${available_camccds[@]}"; do
    cam=${camccd%%:*}
    ccd=${camccd##*:}
    local_camccd="sector${sector}/cam${cam}_ccd${ccd}"
    orbit_dir="$dtarget_base/cam${cam}-ccd${ccd}/$orbit"

    if [[ ! -s "$local_camccd/phot.data" ]]; then
        echo "  WARNING: missing/empty $local_camccd/phot.data; skipping cam${cam}_ccd${ccd}"
        continue
    fi

    echo "  Preparing cam${cam}_ccd${ccd} $orbit"
    cp "$local_camccd/phot.data" "$orbit_dir/phot.data"

    # Clear only this orbit's per-slice temporary LC products before rerunning it.
    for slice_dir in "$orbit_dir"/slice*; do
        [[ -d "$slice_dir" ]] || continue
        rm -rf "$slice_dir/$lcdir" "$slice_dir/bkg_phot/$lcdir"
    done

    do_job=$(sbatch \
        --parsable \
        --array=0-5%6 \
        --time=1-00:00:00 \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        "$script_dir/do_phot_em2.sbatch" \
        "$sector" "$cam" "$ccd" "$orbit" "$lcdir")
    echo "    do_phot cam${cam}_ccd${ccd}: $do_job"

    copy_job=$(sbatch \
        --parsable \
        --array=0-5%6 \
        --time=1-00:00:00 \
        --output="${logdir}/%x.o%A_%a" \
        --error="${logdir}/%x.e%A_%a" \
        --dependency=afterany:"$do_job" \
        "$script_dir/copy_phot_em2.sbatch" \
        "$sector" "$cam" "$ccd" "$orbit" "$lcdir")
    echo "    copy_phot cam${cam}_ccd${ccd}: $copy_job"

    do_jobs+=("$do_job")
    copy_jobs+=("$copy_job")
done

if [[ ${#copy_jobs[@]} -eq 0 ]]; then
    echo "ERROR: no photometry jobs were submitted"
    exit 1
fi
wait_for_jobs "copy_phot" "${copy_jobs[@]}"

echo ""
echo "[4/5] Sorting/deduplicating aggregate raw LC files and purging $orbit slice LC products..."
for camccd in "${available_camccds[@]}"; do
    cam=${camccd%%:*}
    ccd=${camccd##*:}
    outbase="sector${sector}/cam${cam}_ccd${ccd}"
    [[ -d "$outbase/$lcdir" ]] || continue

    find "$outbase/$lcdir" -maxdepth 1 -type f -name 'lc_*.txt' ! -name '*_cleaned' ! -name '*.tmp' -print0 \
        | while IFS= read -r -d '' lcfile; do
            sort -nuk1 "$lcfile" > "${lcfile}.tmp" && mv "${lcfile}.tmp" "$lcfile"
            bkgfile="$outbase/bkg_phot/$lcdir/$(basename "$lcfile")"
            if [[ -f "$bkgfile" ]]; then
                sort -nuk1 "$bkgfile" > "${bkgfile}.tmp" && mv "${bkgfile}.tmp" "$bkgfile"
            fi
        done

    orbit_dir="$dtarget_base/cam${cam}-ccd${ccd}/$orbit"
    for slice_dir in "$orbit_dir"/slice*; do
        [[ -d "$slice_dir" ]] || continue
        rm -rf "$slice_dir/$lcdir" "$slice_dir/bkg_phot/$lcdir"
    done
done

echo ""
echo "[5/5] Running cumulative analysis/plot generation..."
analysis_job=$(sbatch --parsable "$script_dir/analysis.sbatch" "$sector")
echo "  analysis job ID: $analysis_job"
wait_for_jobs "analysis" "$analysis_job"

echo ""
echo "======================================================"
echo " Incremental orbit pipeline finished."
echo " Sector: $sectoruse"
echo " Orbit : $orbit"
echo " Analysis job: $analysis_job"
echo "======================================================"
