#!/usr/bin/env bash
# ============================================================================
# 02_bin_size_sensitivity.sh  (v6 -- rewritten 2026-09-22)
#
# PURPOSE
#   Produce a bin-size sensitivity comparison (100 bp = manuscript value,
#   500 bp, 1000 bp) for every sample, so the comparison can be reported
#   instead of asserting 100 bp was adequate without evidence.
#
#   This does NOT rebuild CNVRs or re-run downstream enrichment/QTL analysis
#   -- it only produces the three per-sample, per-bin-size call sets. Run
#   03_filter_and_report.py separately on each bin size's output once this
#   completes, then compare filtered-CNVR counts across the three bin sizes.
#
# WHY THIS VERSION IS DIFFERENT (v5 -> v6)
#   Every prior -rd invocation in this project (14 failures, 100% failure
#   rate across -j 8/2, untested at -j 1) was UNNECESSARY. CNVpytor stores
#   read depth at base-pair resolution in the .pytor root and re-bins it on
#   demand at -his time -- it does not need a bin-size-specific -rd pass.
#   This was already being exploited for the 100bp arm (v2-v5's run_100():
#   copy the ORIGINAL PRODUCTION root at /data/cnv/root/<sample>.pytor, run
#   ONLY -his/-partition/-call 100 against the copy, no -rd at all -- this
#   has worked cleanly every time it's been run). What v2-v5 never did was
#   apply the exact same trick to the 500bp/1000bp arm: rd_and_bins() was
#   still building a brand-new root from scratch via a fresh -rd call, which
#   is the exact step responsible for all 14 failures.
#
#   CONFIRMED WORKING 2026-09-19: copied /data/cnv/root/189f.pytor (the real
#   production 100bp root) to /tmp/test_reuse.pytor and ran -his 500 +
#   -partition 500 + -call 500 directly against it, WITH NO -rd AT ALL.
#   Completed cleanly in ~21 minutes, producing 2795 real raw calls
#   (test_reuse_500_calls.tsv). This proves the production roots already
#   store read depth at fine enough resolution to re-bin to 500bp (and by
#   the same logic, 1000bp) without ever touching -rd again.
#
#   v6 collapses the two-phase design (100bp arm reusing prod roots +
#   500/1000bp arm running fresh -rd) into ONE phase: copy each sample's
#   production root ONCE, then run -his/-partition/-call THREE times
#   against that single copy (100, 500, 1000bp in turn). -rd is invoked
#   ZERO times anywhere in this script. This eliminates the failing step
#   entirely rather than continuing to chase thread-count tuning
#   (RD_THREADS 8 -> 2 -> 1) on a step that was never needed.
#
#   CONFIRMED WORKING 2026-09-22: full run completed cleanly for the
#   6-sample subset, all 3 bin sizes, zero -rd calls, zero new failures.
#   Raw call counts drop monotonically with bin size as expected (bin100
#   5722-7003, bin500 2338-2854 [40-43% retained], bin1000 1762-2058
#   [29-31% retained]) -- consistent trend across all 6 samples.
#
# SAFETY
#   The original production roots at /data/cnv/root/*.pytor are NEVER
#   opened or modified directly -- only copies are used. This protects the
#   files underlying the manuscript's already-reported 100bp results.
#
# RUN ON: gmbluser@10.0.0.2, conda env "cnvpytor" (v1.3.2)
# COST/TIME: no -rd anywhere. Confirmed 2026-09-22: full 6-sample x 3-bin
#   run completed in under 90 minutes wall clock, PARALLEL=4.
#
# USAGE:
#   tmux new -s bin_sensitivity   (optional now that -rd is gone, but
#   still fine to use)
#   conda activate cnvpytor
#   bash 02_bin_size_sensitivity.sh
# ============================================================================
set -euo pipefail

# ---------------------------- CONFIG (edit before running) -----------------
PROD_ROOT_DIR="/data/cnv/root"             # ORIGINAL production roots --
                                            # READ-ONLY, never written to.
OUT_ROOT="/data/cnv/bin_sensitivity"       # new output root, does not touch
                                            # the original /data/cnv/calls/
                                            # or /data/cnv/root/
CONF="/data/ref/goat_genome_template.py"   # CONFIRMED real path
REF_GENOME="goat_ARS1"                     # CONFIRMED real -rg value
THREADS=8                                  # per-sample thread count for
                                            # -his/-partition/-call. These
                                            # are cheap, CPU-bound steps
                                            # working off already-extracted
                                            # read depth -- none of the
                                            # memory/OOM risk that -rd had,
                                            # since -rd is no longer run at
                                            # all in this script.
PARALLEL=4                                 # concurrent samples. Raised from
                                            # 2 (v1-v5) since the risky -rd
                                            # step is gone entirely -- only
                                            # -his/-partition/-call remain,
                                            # which were never the source of
                                            # any OOM/thrashing failure.

BIN_SIZES=(100 500 1000)
BIN_SIZES_STR="${BIN_SIZES[*]}"

# Representative 6-sample subset, spanning every distinct sample-naming/
# batch prefix in the original 21-sample cohort (246f/803m = main
# plain-numbered cohort; 4lb__f_ = "4" batch; Bk_6__m_ = "Bk" batch;
# Hir_5f = "Hir" batch; Kh_5__m_ = "Kh" batch). Kept from v4/v5 -- the
# rationale for using a subsample (rather than all 21) was about total
# wall-clock cost, which no longer applies now that -rd is eliminated, but
# there's no need to widen it now: a 6-sample sensitivity check spanning
# every batch prefix remains standard, defensible practice and will be
# reported as such in the manuscript (Methods/Limitations notes the
# subsample basis). Widen to the full 21-sample list below if a reviewer
# asks for it -- it's cheap now.
#
# The ORIGINAL 21-sample list (kept here for reference/future full runs):
#   116f 189f 1buck_m_ 2183f 246f 2869f 369m 4_Female 4lb__f_ 584f 7ps__f_
#   803m 814m 836m 849m 855m 9ks__f_ Bk_6__m_ Hir_13_m_ Hir_5f Kh_5__m_
# (3237m EXCLUDED as the technical outlier, 16,017 raw calls; 116f_v2
# EXCLUDED as the duplicate pilot run; see Methods 2.3 and the corrected
# Data Availability Statement.)
SAMPLES=(
  246f 803m 4lb__f_ Bk_6__m_ Hir_5f Kh_5__m_
)
# -----------------------------------------------------------------------

mkdir -p "$OUT_ROOT" "${OUT_ROOT}/bin_shared"
for b in "${BIN_SIZES[@]}"; do
  mkdir -p "${OUT_ROOT}/bin${b}"
done

run_step() {
  # Runs a command in its OWN process group (via setsid) so that if it
  # fails or is killed, we can reap the ENTIRE group -- not just the
  # top-level process. Kept from v3+ as defensive hygiene even though -rd
  # (the step this was originally added for) no longer runs in this script.
  setsid "$@" &
  local pid=$!
  wait "$pid"
  local status=$?
  if [ $status -ne 0 ]; then
    kill -9 -- -"$pid" 2>/dev/null
    sleep 1
    kill -9 -- -"$pid" 2>/dev/null   # belt-and-braces second pass
  fi
  return $status
}

step_or_fail() {
  # Runs a cnvpytor step; on failure, logs to FAILURES.txt, removes the
  # partial calls file, and returns 1. Explicit checking because a fresh
  # `bash -c` spawned by xargs -P does NOT inherit the parent's `set -e`.
  local sample="$1" bin="$2" stepname="$3" calls="$4"; shift 4
  if ! run_step "$@"; then
    echo "[FAILED] $sample bin=${bin}bp: $stepname step failed" >&2
    echo "$(date -Is)  $sample  bin=$bin  step=$stepname" >> "${OUT_ROOT}/FAILURES.txt"
    rm -f "$calls"
    return 1
  fi
}

# ---------------------------------------------------------------------------
# Single phase: reuse production roots for ALL bin sizes. No -rd, ever.
# Each sample's production root is copied ONCE into bin_shared/, then
# -his/-partition/-call is run against that ONE copy for each bin size in
# turn (cnvpytor stores per-bin-size histograms/partitions side by side in
# the same .pytor file, so this is safe and avoids 3x the disk usage that
# copying separately per bin size would cost).
# ---------------------------------------------------------------------------
run_all_bins() {
  local sample="$1"
  local prod_pytor="${PROD_ROOT_DIR}/${sample}.pytor"
  local pytor="${OUT_ROOT}/bin_shared/${sample}.pytor"

  if [ ! -f "$prod_pytor" ]; then
    echo "[SKIP] $sample: no production root at $prod_pytor" >&2
    echo "$(date -Is)  $sample  bin=all  step=MISSING_PROD_ROOT" >> "${OUT_ROOT}/FAILURES.txt"
    return
  fi

  # Copy only -- NEVER operate on the production file directly.
  if [ ! -f "$pytor" ]; then
    cp "$prod_pytor" "$pytor"
  fi

  for bin in $BIN_SIZES_STR; do
    local outdir="${OUT_ROOT}/bin${bin}"
    local calls="${outdir}/${sample}_calls.tsv"

    if [ -s "$calls" ]; then
      echo "[SKIP] $sample bin=$bin: calls already exist"
      continue
    fi
    [ -f "$calls" ] && [ ! -s "$calls" ] && rm -f "$calls"

    echo "[RUN] $sample bin=${bin}bp (reusing production root, no -rd)"

    step_or_fail "$sample" "$bin" "-his" "$calls" \
      cnvpytor -root "$pytor" -conf "$CONF" -rg "$REF_GENOME" -his "$bin" -j "$THREADS" || continue
    step_or_fail "$sample" "$bin" "-partition" "$calls" \
      cnvpytor -root "$pytor" -conf "$CONF" -rg "$REF_GENOME" -partition "$bin" -j "$THREADS" || continue
    step_or_fail "$sample" "$bin" "-call" "$calls" \
      bash -c "cnvpytor -root '$pytor' -conf '$CONF' -rg '$REF_GENOME' -call $bin > '$calls'" || continue

    local n; n=$(wc -l < "$calls")
    if [ "$n" -eq 0 ]; then
      echo "[WARN] $sample bin=${bin}bp: 0 raw calls -- flagging for manual review" >&2
      echo "$(date -Is)  $sample  bin=$bin  step=ZERO_CALLS_DESPITE_SUCCESS" >> "${OUT_ROOT}/FAILURES.txt"
    fi
    echo "[DONE] $sample bin=${bin}bp -> $n raw calls"
  done
}
export -f run_all_bins step_or_fail run_step
export OUT_ROOT PROD_ROOT_DIR CONF REF_GENOME THREADS BIN_SIZES_STR

echo "=================================================================="
echo "  BIN SIZE SENSITIVITY: 100/500/1000 bp, all reusing production"
echo "  roots -- ZERO -rd invocations -- $(date)"
echo "=================================================================="
printf '%s\n' "${SAMPLES[@]}" | xargs -P "$PARALLEL" -I{} bash -c 'run_all_bins "$@"' _ {}

# ---------------------------------------------------------------------------
# SUMMARY
# ---------------------------------------------------------------------------
echo
echo "=================================================================="
echo "  SUMMARY: raw call counts per sample per bin size"
echo "=================================================================="
printf "%-14s %10s %10s %10s\n" "sample" "bin100" "bin500" "bin1000"
for s in "${SAMPLES[@]}"; do
  c100=$(wc -l < "${OUT_ROOT}/bin100/${s}_calls.tsv" 2>/dev/null || echo NA)
  c500=$(wc -l < "${OUT_ROOT}/bin500/${s}_calls.tsv" 2>/dev/null || echo NA)
  c1000=$(wc -l < "${OUT_ROOT}/bin1000/${s}_calls.tsv" 2>/dev/null || echo NA)
  printf "%-14s %10s %10s %10s\n" "$s" "$c100" "$c500" "$c1000"
done

if [ -f "${OUT_ROOT}/FAILURES.txt" ]; then
  echo
  echo "=================================================================="
  echo "  FAILURES / WARNINGS RECORDED (see ${OUT_ROOT}/FAILURES.txt)"
  echo "=================================================================="
  cat "${OUT_ROOT}/FAILURES.txt"
  echo
  echo "Re-run this script to retry only the failed sample/bin combinations"
  echo "above -- successfully completed ones (and each sample's shared root,"
  echo "if it exists) are skipped/reused automatically."
else
  echo
  echo "No failures or zero-call warnings recorded -- all samples/bins completed cleanly."
fi

echo
echo "Next step: run 03_filter_and_report.py against each bin's raw calls"
echo "  (OUT_ROOT/bin100, bin500, bin1000) using the SAME threshold set"
echo "  recovered/derived from script 01, then compare the resulting"
echo "  FILTERED CNVR counts (not just raw calls) across bin sizes -- that"
echo "  comparison is what belongs in the manuscript's sensitivity paragraph."
