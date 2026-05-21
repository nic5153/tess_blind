#!/usr/bin/env bash
#SBATCH --job-name="agg_errors"
#SBATCH --output=/lustre/research/mfausnau/logs/%x.o%A
#SBATCH --error=/lustre/research/mfausnau/logs/%x.e%A
#SBATCH --partition=nocona
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-0:30:00

set -eo pipefail

logdir="${LOG_DIR:-/lustre/research/mfausnau/logs}"
dhome=$(pwd)

# Arguments: <sector> <cam> <ccd> <lcdir> <do_phot_jobid1> [...] -- <copy_phot_jobid1> [...] -- <clean_phot_jobid1> [...]
# For submit_all: cam and ccd will be "all"
if [[ $# -lt 5 ]]; then
    echo "Usage: $0 <sector> <cam> <ccd> <lcdir> <do_phot_jobid1> [...] -- <copy_phot_jobid1> [...] -- <clean_phot_jobid1> [...]"
    exit 1
fi

sector=$1
cam=$2
ccd=$3
lcdir=$4
shift 4

# Split remaining args on "--" into do_phot, copy_phot, and clean_phot job ID lists
do_phot_jobids=()
copy_phot_jobids=()
clean_phot_jobids=()
target=do_phot_jobids
for arg in "$@"; do
    if [[ "$arg" == "--" ]]; then
        if [[ "$target" == "do_phot_jobids" ]]; then
            target=copy_phot_jobids
        else
            target=clean_phot_jobids
        fi
    else
        eval "${target}+=(\"\$arg\")"
    fi
done

sectoruse=$(printf "s%04d" "$sector")
timestamp=$(date +%Y%m%d_%H%M%S)
outfile="$dhome/errorlog-${sectoruse}-cam${cam}-ccd${ccd}-${timestamp}.txt"

echo "Aggregating errors for sector=$sectoruse cam=$cam ccd=$ccd"
echo "do_phot job IDs:    ${do_phot_jobids[*]}"
echo "copy_phot job IDs:  ${copy_phot_jobids[*]}"
echo "clean_phot job IDs: ${clean_phot_jobids[*]}"
echo "Output: $outfile"

{
    echo "========================================================"
    echo " Error report for sector=$sectoruse cam=$cam ccd=$ccd"
    echo " lcdir=$lcdir"
    echo " Generated: $(date)"
    echo " do_phot job IDs:    ${do_phot_jobids[*]}"
    echo " copy_phot job IDs:  ${copy_phot_jobids[*]}"
    echo " clean_phot job IDs: ${clean_phot_jobids[*]}"
    echo "========================================================"
    echo ""

    total_failed=0
    do_errors_found=0

    echo "--- do_phot failures ---"
    echo ""

    for jobid in "${do_phot_jobids[@]}"; do
        mapfile -t efiles < <(ls "$logdir/do_phot.e${jobid}_"* 2>/dev/null || true)

        if [[ ${#efiles[@]} -eq 0 ]]; then
            echo "WARNING: No error log files found for do_phot job $jobid"
            continue
        fi

        for efile in "${efiles[@]}"; do
            [[ -s "$efile" ]] || continue

            task_id=$(basename "$efile" | grep -oP '(?<=_)\d+$')
            echo "--- do_phot job ${jobid}, task ${task_id} ---"
            cat "$efile"
            echo ""
            (( do_errors_found++ )) || true

            # For each phot2.csh failure, also print the dates file
            while IFS= read -r line; do
                slice_path=$(echo "$line" | grep -oP '(?<=in ).*(?= \(exit)')
                [[ -z "$slice_path" ]] && continue
                base_slice="${slice_path%/bkg_phot}"
                dates_file="$base_slice/dates"
                if [[ -e "$dates_file" ]]; then
                    echo "  Dates for $base_slice:"
                    sed 's/^/    /' "$dates_file"
                else
                    echo "  (dates file not found in $base_slice)"
                fi
                echo ""
                (( total_failed++ )) || true
            done < <(grep "^FAILED" "$efile" 2>/dev/null || true)
        done
    done

    echo "--- copy_phot errors ---"
    echo ""

    # ---- copy_phot errors (content of error logs only, no dates) ----
    copy_errors_found=0
    for jobid in "${copy_phot_jobids[@]}"; do
        mapfile -t efiles < <(ls "$logdir/copy_phot.e${jobid}_"* 2>/dev/null || true)

        if [[ ${#efiles[@]} -eq 0 ]]; then
            echo "WARNING: No error log files found for copy_phot job $jobid"
            continue
        fi

        for efile in "${efiles[@]}"; do
            # Skip empty error logs
            if [[ ! -s "$efile" ]]; then
                continue
            fi
            task_id=$(basename "$efile" | grep -oP '(?<=_)\d+$')
            echo "--- copy_phot job ${jobid}, task ${task_id} ---"
            cat "$efile"
            echo ""
            (( copy_errors_found++ )) || true
        done
    done

    if [[ $copy_errors_found -eq 0 ]]; then
        echo "(No copy_phot errors found)"
        echo ""
    fi

    echo "--- clean_phot errors ---"
    echo ""

    clean_errors_found=0
    for jobid in "${clean_phot_jobids[@]}"; do
        mapfile -t efiles < <(ls "$logdir/clean_phot.e${jobid}_"* 2>/dev/null || true)

        if [[ ${#efiles[@]} -eq 0 ]]; then
            echo "WARNING: No error log files found for clean_phot job $jobid"
            continue
        fi

        for efile in "${efiles[@]}"; do
            [[ -s "$efile" ]] || continue
            task_id=$(basename "$efile" | grep -oP '(?<=_)\d+$')
            echo "--- clean_phot job ${jobid}, task ${task_id} ---"
            cat "$efile"
            echo ""
            (( clean_errors_found++ )) || true
        done
    done

    if [[ $clean_errors_found -eq 0 ]]; then
        echo "(No clean_phot errors found)"
        echo ""
    fi

    echo "========================================================"
    echo " do_phot tasks with errors:    $do_errors_found"
    echo " Total failed phot2.csh calls: $total_failed"
    echo " copy_phot tasks with errors:  $copy_errors_found"
    echo " clean_phot tasks with errors: $clean_errors_found"
    echo "========================================================"

} > "$outfile"

echo "Done. $total_failed failure(s) reported."
cat "$outfile"
