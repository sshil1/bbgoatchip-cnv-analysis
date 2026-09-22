#!/usr/bin/env bash
# ============================================================================
# 01_recover_historical_thresholds.sh
#
# PURPOSE
#   Before re-deriving CNVpytor filter thresholds from scratch, check whether
#   the ORIGINAL filtering run (the one that already produced
#   population_cnvrs_filtered_clean_nuclear.bed, 9,717 CNVRs) left any trace
#   of the actual e-val1 / q0 / pN / dG cutoffs used. If it did, those are the
#   numbers that belong in Methods 2.2 and the supplementary QC table --
#   re-deriving a DIFFERENT threshold set now would silently change the
#   9,717-CNVR result already reported throughout the manuscript.
#
#   Run this FIRST. Only fall back to 03_filter_and_report.py's default
#   threshold set if this script finds nothing.
#
# RUN ON: gmbluser@10.0.0.2 (has internet, original CNVpytor pipeline)
# USAGE:
#   tmux new -s recover_thresholds
#   bash 01_recover_historical_thresholds.sh 2>&1 | tee recover_thresholds_report.txt
#   # detach: Ctrl-b d ; reattach: tmux attach -t recover_thresholds
# ============================================================================
set -uo pipefail

echo "=== 1. Shell history (all users this account can read) ==="
for h in ~/.bash_history ~/.zsh_history; do
  [ -f "$h" ] && grep -nE "cnvpytor|eval1|q0|pN|dG" "$h" | grep -iE "call|filter"
done

echo
echo "=== 2. Any script files under home dir or /data referencing CNVpytor thresholds ==="
find "$HOME" /data /opt -maxdepth 6 \( -iname "*.py" -o -iname "*.sh" -o -iname "*.R" \) \
     -newermt "2026-06-01" 2>/dev/null \
  | xargs grep -lE "eval1|e_val1|e-val1|pN|dG|q0" 2>/dev/null

echo
echo "=== 3. Grep those hits for the actual numeric comparisons (e.g. eval1 < 0.05) ==="
find "$HOME" /data /opt -maxdepth 6 \( -iname "*.py" -o -iname "*.sh" -o -iname "*.R" \) \
     -newermt "2026-06-01" 2>/dev/null \
  | xargs grep -nE "(eval1|e_val1|e-val1|pN|dG|q0)\s*[<>=]" 2>/dev/null

echo
echo "=== 4. nohup / log files near the CNV calls directory ==="
find /data/cnv -iname "*.log" -o -iname "nohup.out" -o -iname "*.out" 2>/dev/null \
  | while read -r f; do
      echo "--- $f ---"
      grep -nE "cnvpytor|eval1|filter|threshold" "$f" 2>/dev/null | head -20
    done

echo
echo "=== 5. Any README / notes files documenting the pipeline ==="
find "$HOME" /data/cnv -maxdepth 4 -iname "README*" -o -iname "*NOTES*" -o -iname "*methods*" 2>/dev/null

echo
echo "=== 6. Compare candidate raw call counts against the known target numbers ==="
echo "    (5,682-7,003 calls/sample except 3237m=16,017; final = 9,717 filtered CNVRs, 21 samples)"
echo "    If a script's output count matches these, its thresholds are almost certainly the real ones."
find /data -iname "*calls*.tsv" -o -iname "*calls*.txt" 2>/dev/null | head -30

echo
echo "=== DONE ==="
echo "If sections 1-5 recovered nothing, the thresholds were most likely applied"
echo "interactively (not scripted) or CNVpytor's tool defaults were used unmodified."
echo "In that case run 03_filter_and_report.py with --thresholds default and clearly"
echo "label the Methods/Supplementary table as 'CNVpytor default quality filters',"
echo "not as breed/study-specific tuned values."
