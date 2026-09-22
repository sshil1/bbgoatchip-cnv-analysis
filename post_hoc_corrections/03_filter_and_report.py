#!/usr/bin/env python3
"""
03_filter_and_report.py

PURPOSE
  CONFIRMED (2026-09-18, via 01_recover_historical_thresholds.sh + a direct
  `wc -l` check against /data/cnv/calls/*.pytor.txt): no per-call e-val1/
  q0/pN/dG quality filter was ever applied to the original CNVpytor output.
  The raw per-sample call counts (5,682-7,003, except 3237m's 16,017) match
  the numbers already reported in the manuscript's QC section exactly, so
  those ARE the unfiltered calls -- they were carried straight into CNVR
  construction, and "filtered CNVRs" in the manuscript refers only to the
  population-level proximity-merge + singleton-removal step, not a call
  quality filter. Methods 2.2 and Limitations have already been corrected
  in build.js to state this accurately instead of the old placeholder.

  This script is therefore now OPTIONAL, not a blocking gap: it applies a
  post hoc per-call quality filter (community-standard CNVpytor thresholds
  by default) as a step-by-step funnel, so you can see what a stricter,
  quality-filtered CNVR set would look like as a strengthening follow-up
  analysis -- e.g. to check whether the headline findings (BMPR1B, the
  Chr.1 NK-cytotoxicity cluster) survive filtering, or to add a genuinely
  new "quality-filtered sensitivity" result to the manuscript.

  CNVpytor -call output columns (no header row in the raw file):
    1  CNV_type        "deletion" | "duplication"
    2  coordinates      chr:start-end
    3  size             bp
    4  cnv_level        normalized read depth in the region (~0.5 = het del,
                         ~0 = hom del, ~1.5 = het dup, ~2 = hom dup, on a
                         diploid-normalized scale)
    5  e_val1           t-test p-value, region vs global read depth
    6  e_val2           Gaussian-likelihood p-value
    7  e_val3           e_val1 computed over the whole (not just call) region
    8  e_val4           e_val2 computed over the whole region
    9  q0               fraction of reads in region with mapping quality 0
    10 pN               fraction of N (unmappable) bases in the region
    11 dG               distance (bp) from region boundary to nearest
                         assembly gap

RUN ON: gmbluser@10.0.0.2
USAGE (confirmed real path on this server):
  python3 03_filter_and_report.py \
      --calls-dir /data/cnv/calls \
      --out-dir   /data/cnv/calls/filter_report \
      --eval1 0.05 --q0 0.5 --pn 0.5 --dg 100000

  NOTE: --calls-dir points straight at /data/cnv/calls -- that directory
  contains the 23 raw *.pytor.txt files (21 analysis samples + the excluded
  3237m and 116f_v2) alongside nothing else with a .pytor.txt extension, so
  no need to copy them elsewhere first.

  The --eval1/--q0/--pn/--dg defaults are standard CNVpytor community
  filtering thresholds (not values recovered from this study's history --
  none exist, per the confirmation above). Report them explicitly as a
  NEW post hoc filter you are applying now, not as thresholds used in the
  original analysis.
"""
import argparse
import csv
import glob
import os
import sys
from collections import OrderedDict

COLS = [
    "cnv_type", "coordinates", "size", "cnv_level",
    "e_val1", "e_val2", "e_val3", "e_val4", "q0", "pN", "dG",
]


def load_calls(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) < 11:
                continue
            row = dict(zip(COLS, parts))
            rows.append(row)
    return rows


def to_float(x, default=None):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def funnel_one_sample(rows, eval1_max, q0_max, pn_max, dg_min):
    """Apply filters sequentially, recording the survivor count after each
    step (order matches the order named in Methods 2.2: e-val1, q0, pN, dG).
    """
    steps = OrderedDict()
    steps["raw"] = rows

    step1 = [r for r in rows if to_float(r["e_val1"], 1.0) < eval1_max]
    steps[f"e_val1<{eval1_max}"] = step1

    step2 = [r for r in step1 if to_float(r["q0"], 1.0) < q0_max]
    steps[f"q0<{q0_max}"] = step2

    step3 = [r for r in step2 if to_float(r["pN"], 1.0) < pn_max]
    steps[f"pN<{pn_max}"] = step3

    step4 = [r for r in step3 if to_float(r["dG"], 0.0) > dg_min]
    steps[f"dG>{dg_min}"] = step4

    return steps


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--calls-dir", required=True,
                     help="Directory of raw per-sample *.pytor.txt "
                          "(CNVpytor -call output, no header). Confirmed "
                          "real path/naming on gmbluserserver: "
                          "/data/cnv/calls/<sample>.pytor.txt")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--eval1", type=float, default=0.05)
    ap.add_argument("--q0", type=float, default=0.5)
    ap.add_argument("--pn", type=float, default=0.5)
    ap.add_argument("--dg", type=float, default=100000)
    ap.add_argument("--pattern", default="*.pytor.txt")
    ap.add_argument("--suffix", default=".pytor.txt",
                     help="Suffix stripped from each matched filename to "
                          "recover the sample name (default matches "
                          "<sample>.pytor.txt).")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    files = sorted(glob.glob(os.path.join(args.calls_dir, args.pattern)))
    if not files:
        sys.exit(f"No files matching {args.pattern} in {args.calls_dir}")

    per_sample_summary = []
    step_names = None

    for f in files:
        sample = os.path.basename(f)
        if sample.endswith(args.suffix):
            sample = sample[: -len(args.suffix)]
        rows = load_calls(f)
        steps = funnel_one_sample(rows, args.eval1, args.q0, args.pn, args.dg)
        if step_names is None:
            step_names = list(steps.keys())

        counts = {k: len(v) for k, v in steps.items()}
        counts["sample"] = sample
        per_sample_summary.append(counts)

        # write the final filtered calls for this sample (for downstream
        # CNVR re-construction if you re-run the merge step with the
        # now-confirmed thresholds)
        final_rows = steps[step_names[-1]]
        with open(os.path.join(args.out_dir, f"{sample}_filtered.tsv"), "w") as out:
            w = csv.writer(out, delimiter="\t")
            for r in final_rows:
                w.writerow([r[c] for c in COLS])

    # ---- per-sample funnel table ----
    funnel_path = os.path.join(args.out_dir, "per_sample_filter_funnel.tsv")
    with open(funnel_path, "w") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["sample"] + step_names)
        for row in per_sample_summary:
            w.writerow([row["sample"]] + [row[s] for s in step_names])

    # ---- pooled totals (this is the number set that goes into the
    # manuscript's Methods 2.2 / supplementary QC table) ----
    totals = {s: sum(r[s] for r in per_sample_summary) for s in step_names}
    totals_path = os.path.join(args.out_dir, "pooled_filter_funnel.tsv")
    with open(totals_path, "w") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["step", "total_calls_all_samples", "n_removed_this_step"])
        prev = None
        for s in step_names:
            removed = "" if prev is None else (totals[prev] - totals[s])
            w.writerow([s, totals[s], removed])
            prev = s

    print(f"Thresholds used: e_val1<{args.eval1}, q0<{args.q0}, "
          f"pN<{args.pn}, dG>{args.dg}")
    print(f"Per-sample funnel written to: {funnel_path}")
    print(f"Pooled funnel written to:     {totals_path}")
    print(f"Filtered per-sample calls in: {args.out_dir}/*_filtered.tsv")
    print()
    print("Paste the pooled_filter_funnel.tsv contents directly into a")
    print("supplementary QC table, and quote the four threshold values")
    print("above verbatim in Methods 2.2 in place of the current")
    print('"being compiled from pipeline run logs" placeholder.')


if __name__ == "__main__":
    main()
