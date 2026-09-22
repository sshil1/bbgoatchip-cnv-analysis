#!/bin/bash
# Runs run_one_sample.sh across every dedup'd goat sample, 2 at a time.
# Usage: run inside a tmux window, in the `cnvpytor` conda env.

set -u

BAM_DIR="/data/bam/dedup"
CALLS_DIR="/data/cnv/calls"
LOG_DIR="/data/logs/cnvpytor"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$CALLS_DIR" "$LOG_DIR"

# Build sample list from every dedup BAM present
SAMPLES=()
for f in "${BAM_DIR}"/*.dedup.bam; do
  SAMPLES+=("$(basename "$f" .dedup.bam)")
done

echo "Found ${#SAMPLES[@]} dedup'd samples:"
printf '  %s\n' "${SAMPLES[@]}"
echo ""
echo "Starting batch (2 samples in parallel, ~30 chromosomes each, ~3h/sample) ..."

printf '%s\n' "${SAMPLES[@]}" | xargs -I{} -P 2 "${SCRIPT_DIR}/run_one_sample.sh" {}

echo ""
echo "=== BATCH COMPLETE $(date) ==="
echo "Summary:"
for s in "${SAMPLES[@]}"; do
  c="${CALLS_DIR}/${s}.pytor.txt"
  if [[ -s "$c" ]]; then
    echo "  $s: $(wc -l < "$c") calls"
  else
    echo "  $s: NO OUTPUT (check ${LOG_DIR}/${s}.log)"
  fi
done
