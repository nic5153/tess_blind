#!/usr/bin/env bash

#SBATCH --job-name=full_tess_pipeline
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --time=3-00:00:00

set -eo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch $0 <sector> [lcdir]"
    echo "Example: sbatch $0 101 lc"
    exit 1
fi

sector=$1
lcdir=${2:-lc}
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source /home/${USER}/setup_vars.sh
logdir="${LOG_DIR:-/lustre/work/nimcclur/TESS/photometry/logs}"
mkdir -p "$logdir"

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

        if [[ "$running" -eq 0 ]]; then
            break
        fi

        sleep 60
    done
    echo "Finished waiting for $label"
}

echo "======================================================"
echo " Full TESS blind pipeline"
echo " Sector : $sector"
echo " lcdir  : $lcdir"
echo " CWD    : $(pwd)"
echo " Scripts: $script_dir"
echo " Logs   : $logdir"
echo "======================================================"

echo ""
echo "[1/6] Creating local sector/cam/ccd/orbit directories..."
mkdir -p "sector${sector}"
for cam in {1..4}; do
    for ccd in {1..4}; do
        for orbit in o1a o1b o2a o2b; do
            mkdir -p "sector${sector}/cam${cam}_ccd${ccd}/${orbit}"
        done
    done
done

echo "  Copying source-stage helper scripts into sector${sector}/"
cp "$script_dir/do_photutils_extract.py" "sector${sector}/do_photutils_extract.py"
cp "$script_dir/rms_source.py" "sector${sector}/rms_source.py"

echo ""
echo "[2/6] Submitting daily RMS generation for all cam/ccd combinations..."
rms_jobs=()
for cam in {1..4}; do
    for ccd in {1..4}; do
        rms_job=$(sbatch --parsable "$script_dir/make_rms_sbatch.sh" "$sector" "$cam" "$ccd")
        echo "  daily_rms cam${cam}_ccd${ccd}: $rms_job"
        rms_jobs+=("$rms_job")
    done
done
wait_for_jobs "daily RMS generation" "${rms_jobs[@]}"

echo ""
echo "[3/6] Running source extraction from daily RMS products..."
source_job=$(sbatch --parsable "$script_dir/source.sbatch" "$sector")
echo "  source job ID: $source_job"
wait_for_jobs "source extraction" "$source_job"

echo ""
echo "[4/6] Submitting photometry pipeline..."
submit_log="full_pipeline_submit_all_sector${sector}_$(date +%Y%m%d_%H%M%S).log"
set +e
bash "$script_dir/submit_all.sh" "$sector" "$lcdir" 2>&1 | tee "$submit_log"
submit_status=${PIPESTATUS[0]}
set -e
if [[ "$submit_status" -ne 0 ]]; then
    echo "ERROR: submit_all.sh failed with status $submit_status"
    exit "$submit_status"
fi

purge_job=$(awk '/purge_lc job ID:/ {print $NF}' "$submit_log" | tail -n 1)
if [[ -z "$purge_job" ]]; then
    echo "ERROR: could not parse purge_lc job ID from $submit_log"
    exit 1
fi
echo "  purge_lc terminal job ID: $purge_job"
wait_for_jobs "photometry pipeline through purge_lc" "$purge_job"

echo ""
echo "[5/6] Untarring image products for plotting..."
untar_job=$(sbatch --parsable "$script_dir/images_untar.sbatch" "$sector")
echo "  images_untar job ID: $untar_job"
wait_for_jobs "image untar" "$untar_job"

echo ""
echo "[6/6] Running analysis and plot generation..."
analysis_job=$(sbatch --parsable "$script_dir/analysis.sbatch" "$sector")
echo "  analysis job ID: $analysis_job"
wait_for_jobs "analysis" "$analysis_job"

echo ""
echo "======================================================"
echo " Full pipeline finished submitting and waiting."
echo " Sector: $sector"
echo " Analysis job: $analysis_job"
echo " Check plots under sector${sector}/cam*_ccd*/plots and sector${sector}/{single,multiple,flag}."
echo "======================================================"
