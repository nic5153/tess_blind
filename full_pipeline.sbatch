#!/usr/bin/env bash

#SBATCH --job-name=full_tess_pipeline
#SBATCH --output=%x.o%j
#SBATCH --error=%x.e%j
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
<<<<<<< HEAD
#SBATCH --time=00:08:00
=======
#SBATCH --time=1-00:00:00
>>>>>>> origin/codex/raw-lc-daily-rms

set -eo pipefail

if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: sbatch $0 <sector> [lcdir]"
    echo "Example: sbatch $0 101 lc"
    exit 1
fi

sector=$1
lcdir=${2:-lc}
if [[ -n "${SLURM_SUBMIT_DIR:-}" && -f "${SLURM_SUBMIT_DIR}/do_photutils_extract.py" ]]; then
    script_dir="$SLURM_SUBMIT_DIR"
else
    script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi

source /home/${USER}/setup_vars.sh
logdir="${LOG_DIR:-/lustre/work/nimcclur/TESS/photometry/logs}"
mkdir -p "$logdir"
<<<<<<< HEAD
=======
ANALYSIS_CPUS=${ANALYSIS_CPUS:-4}
RUN_MULTIPLE_EXPORT=${RUN_MULTIPLE_EXPORT:-1}
RUN_CATALOG_ENRICHMENT=${RUN_CATALOG_ENRICHMENT:-0}
CATALOG_SLEEP=${CATALOG_SLEEP:-0.4}
CATALOG_LIMIT=${CATALOG_LIMIT:-}
MULTIPLE_CSV=${MULTIPLE_CSV:-multiple_candidates.csv}
MULTIPLE_HTML=${MULTIPLE_HTML:-multiple_candidates.html}
MULTIPLE_ENRICHED_CSV=${MULTIPLE_ENRICHED_CSV:-multiple_candidates_enriched.csv}
>>>>>>> origin/codex/raw-lc-daily-rms

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
<<<<<<< HEAD
            if squeue -h -j "$job" | grep -q .; then
=======
            if squeue -h -j "$job" 2>/dev/null | grep -q .; then
>>>>>>> origin/codex/raw-lc-daily-rms
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
<<<<<<< HEAD
echo "======================================================"

echo ""
echo "[1/6] Creating local sector/cam/ccd/orbit directories..."
=======
echo " Analysis CPUs : $ANALYSIS_CPUS"
echo " Multiple export: $RUN_MULTIPLE_EXPORT"
echo " Catalog enrich : $RUN_CATALOG_ENRICHMENT"
echo "======================================================"

echo ""
echo "[1/7] Creating local sector/cam/ccd/orbit directories..."
>>>>>>> origin/codex/raw-lc-daily-rms
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
<<<<<<< HEAD
echo "[2/6] Submitting daily RMS generation for all cam/ccd combinations..."
=======
echo "[2/7] Submitting daily RMS generation for all cam/ccd combinations..."
>>>>>>> origin/codex/raw-lc-daily-rms
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
<<<<<<< HEAD
echo "[3/6] Running source extraction from daily RMS products..."
=======
echo "[3/7] Running source extraction from daily RMS products..."
>>>>>>> origin/codex/raw-lc-daily-rms
source_job=$(sbatch --parsable "$script_dir/source.sbatch" "$sector")
echo "  source job ID: $source_job"
wait_for_jobs "source extraction" "$source_job"

echo ""
<<<<<<< HEAD
echo "[4/6] Submitting photometry pipeline..."
=======
echo "[4/7] Submitting photometry pipeline..."
>>>>>>> origin/codex/raw-lc-daily-rms
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

<<<<<<< HEAD
echo ""
echo "[5/6] Untarring image products for plotting..."
=======
lc_count=$(find "sector${sector}" -path "*/${lcdir}/lc_*.txt" -type f -printf '.' | wc -c)
bkg_lc_count=$(find "sector${sector}" -path "*/bkg_phot/${lcdir}/lc_*.txt" -type f -printf '.' | wc -c)
echo "Copied source light curves: $lc_count"
echo "Copied background light curves: $bkg_lc_count"

if [[ "$lc_count" -eq 0 || "$bkg_lc_count" -eq 0 ]]; then
    echo "ERROR: photometry/copy stage produced no local light curves; aborting before image untar/analysis." >&2
    echo "Check ${logdir}/do_phot.e*, ${logdir}/copy_phot.e*, and ${logdir}/clean_phot.e*." >&2
    exit 1
fi

echo ""
echo "[5/7] Untarring image products for plotting..."
>>>>>>> origin/codex/raw-lc-daily-rms
untar_job=$(sbatch --parsable "$script_dir/images_untar.sbatch" "$sector")
echo "  images_untar job ID: $untar_job"
wait_for_jobs "image untar" "$untar_job"

echo ""
<<<<<<< HEAD
echo "[6/6] Running analysis and plot generation..."
analysis_job=$(sbatch --parsable "$script_dir/analysis.sbatch" "$sector")
=======
echo "[6/7] Running analysis and plot generation..."
analysis_job=$(sbatch --parsable --cpus-per-task="$ANALYSIS_CPUS" "$script_dir/analysis.sbatch" "$sector")
>>>>>>> origin/codex/raw-lc-daily-rms
echo "  analysis job ID: $analysis_job"
wait_for_jobs "analysis" "$analysis_job"

echo ""
<<<<<<< HEAD
=======
echo "[7/7] Exporting multiple-candidate CSV/HTML products..."
if [[ "$RUN_MULTIPLE_EXPORT" == "1" ]]; then
    if python "$script_dir/export_multiple_candidates.py" \
        --root "." \
        --outfile "$MULTIPLE_CSV" \
        --html "$MULTIPLE_HTML" \
        --enriched-csv "$MULTIPLE_ENRICHED_CSV"; then
        echo "  Wrote base multiple-candidate products."
    else
        echo "WARNING: export_multiple_candidates.py failed; continuing." >&2
    fi

    if [[ "$RUN_CATALOG_ENRICHMENT" == "1" && -s "$MULTIPLE_CSV" ]]; then
        enrich_args=("$MULTIPLE_CSV" --output "$MULTIPLE_ENRICHED_CSV" --sleep "$CATALOG_SLEEP")
        if [[ -n "$CATALOG_LIMIT" ]]; then
            enrich_args+=(--limit "$CATALOG_LIMIT")
        fi

        if python "$script_dir/enrich_multiple_candidates.py" "${enrich_args[@]}"; then
            python "$script_dir/export_multiple_candidates.py" \
                --root "." \
                --outfile "$MULTIPLE_CSV" \
                --html "$MULTIPLE_HTML" \
                --enriched-csv "$MULTIPLE_ENRICHED_CSV" \
                || echo "WARNING: enriched HTML refresh failed; continuing." >&2
        else
            echo "WARNING: enrich_multiple_candidates.py failed; base CSV/HTML remain available." >&2
        fi
    else
        echo "  Catalog enrichment skipped."
    fi
else
    echo "  Multiple-candidate export skipped."
fi

echo ""
>>>>>>> origin/codex/raw-lc-daily-rms
echo "======================================================"
echo " Full pipeline finished submitting and waiting."
echo " Sector: $sector"
echo " Analysis job: $analysis_job"
<<<<<<< HEAD
=======
echo " Multiple CSV: $MULTIPLE_CSV"
echo " Multiple HTML: $MULTIPLE_HTML"
>>>>>>> origin/codex/raw-lc-daily-rms
echo " Check plots under sector${sector}/cam*_ccd*/plots and sector${sector}/{single,multiple,flag}."
echo "======================================================"
