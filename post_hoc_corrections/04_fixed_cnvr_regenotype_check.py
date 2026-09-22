#!/usr/bin/env python3
"""
04_fixed_cnvr_regenotype_check.py

PURPOSE
  Answers the audit's "fixed CNVR classification unconfirmed by per-sample
  re-genotyping" concern directly.

  As currently constructed, a CNVR is labeled "fixed" (21/21 support) if
  each of the 21 samples has *some* call that was merged into it during the
  <=1kb proximity-merge step. That merge step only requires positional
  proximity between per-sample calls, not that each sample's individual call
  actually reciprocally overlaps most of the final CNVR's footprint. Two
  samples whose individual deletions each cover a different 20% sliver of a
  large CNVR, on opposite ends, would both count toward "21/21 support" under
  the proximity-merge rule while genuinely disagreeing about where (or how
  much of) the deletion actually is.

  This script re-checks every CNVR currently labeled "fixed" using RECIPROCAL
  OVERLAP instead of simple proximity: for each sample, does that sample have
  an individual call whose overlap with the CNVR is >= X% of BOTH the call's
  length AND the CNVR's length (default X=50, i.e. standard reciprocal
  overlap used in SV benchmarking, e.g. Truvari/SVanalyzer conventions)?

  Output: for every currently-"fixed" CNVR, the TRUE reciprocal-overlap
  support count out of 21. Anything that drops below 21 should be relabeled
  in the manuscript (e.g. downgraded from "fixed" to "near-fixed, N/21 by
  reciprocal overlap"), and the 254-gene FIXED-set enrichment analysis
  (NK-cytotoxicity, p=1.5e-21) should be re-run on whatever gene list results
  from the corrected, stricter "fixed" set.

INPUTS
  --population-bed   population_cnvrs_filtered_clean_nuclear.bed
                      Expected columns (tab- or space-delimited, no header):
                      chrom  start  end  cnv_type  n_support  sample_list_csv
                      (matches the format already used in this project's
                      cnvr_qtl_overlaps.bed-style files)
  --per-sample-dir   directory of per-sample FILTERED calls, i.e. the
                      *_filtered.tsv files produced by
                      03_filter_and_report.py (coordinates column format
                      "chr:start-end")

OUTPUT
  fixed_cnvr_regenotype_report.tsv -- one row per CNVR originally labeled
  fixed (n_support == n_total_samples), with:
    chrom, start, end, cnv_type, claimed_support, reciprocal_support,
    reciprocal_pct_of_claimed, still_fixed (bool), samples_failing

USAGE
  python3 04_fixed_cnvr_regenotype_check.py \
      --population-bed /data/cnv/calls/population_cnvrs_filtered_clean_nuclear.bed \
      --per-sample-dir /data/cnv/calls/filter_report \
      --n-total-samples 21 \
      --min-reciprocal-pct 50 \
      --out fixed_cnvr_regenotype_report.tsv
"""
import argparse
import csv
from collections import defaultdict


def parse_population_bed(path):
    cnvrs = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split("\t") if "\t" in line else line.split()
            if len(parts) < 6:
                continue
            chrom, start, end, cnv_type, n_support, sample_list = parts[:6]
            cnvrs.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "cnv_type": cnv_type,
                "n_support": int(n_support),
                "samples": [s for s in sample_list.split(",") if s],
            })
    return cnvrs


def parse_sample_calls(per_sample_dir):
    """Returns dict: sample_name -> list of (chrom, start, end, cnv_type)"""
    import glob
    import os
    calls = defaultdict(list)
    for f in glob.glob(os.path.join(per_sample_dir, "*_filtered.tsv")):
        sample = os.path.basename(f).replace("_filtered.tsv", "")
        with open(f) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                cnv_type, coord = parts[0], parts[1]
                try:
                    chrom, rng = coord.split(":")
                    start, end = rng.split("-")
                    calls[sample].append((chrom, int(start), int(end), cnv_type))
                except ValueError:
                    continue
    return calls


def reciprocal_overlap_ok(a_start, a_end, b_start, b_end, min_pct):
    ov_start = max(a_start, b_start)
    ov_end = min(a_end, b_end)
    if ov_end <= ov_start:
        return False, 0.0
    ov = ov_end - ov_start
    a_len = a_end - a_start
    b_len = b_end - b_start
    pct_a = 100.0 * ov / a_len if a_len > 0 else 0
    pct_b = 100.0 * ov / b_len if b_len > 0 else 0
    ok = pct_a >= min_pct and pct_b >= min_pct
    return ok, min(pct_a, pct_b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--population-bed", required=True)
    ap.add_argument("--per-sample-dir", required=True)
    ap.add_argument("--n-total-samples", type=int, default=21)
    ap.add_argument("--min-reciprocal-pct", type=float, default=50.0)
    ap.add_argument("--out", default="fixed_cnvr_regenotype_report.tsv")
    args = ap.parse_args()

    cnvrs = parse_population_bed(args.population_bed)
    sample_calls = parse_sample_calls(args.per_sample_dir)

    fixed = [c for c in cnvrs if c["n_support"] >= args.n_total_samples]
    print(f"Loaded {len(cnvrs)} total CNVRs; {len(fixed)} currently labeled "
          f"'fixed' ({args.n_total_samples}/{args.n_total_samples} support).")
    print(f"Loaded per-sample call sets for {len(sample_calls)} samples "
          f"from {args.per_sample_dir}")

    rows = []
    n_downgraded = 0
    for c in fixed:
        failing = []
        recip_support = 0
        for sample in c["samples"]:
            calls = sample_calls.get(sample, [])
            best = 0.0
            hit = False
            for (chrom, s, e, ctype) in calls:
                if chrom != c["chrom"] or ctype != c["cnv_type"]:
                    continue
                ok, pct = reciprocal_overlap_ok(c["start"], c["end"], s, e,
                                                 args.min_reciprocal_pct)
                best = max(best, pct)
                if ok:
                    hit = True
                    break
            if hit:
                recip_support += 1
            else:
                failing.append(f"{sample}(best={best:.0f}%)")

        still_fixed = recip_support >= args.n_total_samples
        if not still_fixed:
            n_downgraded += 1
        rows.append({
            "chrom": c["chrom"], "start": c["start"], "end": c["end"],
            "cnv_type": c["cnv_type"], "claimed_support": c["n_support"],
            "reciprocal_support": recip_support,
            "reciprocal_pct_of_claimed": round(
                100.0 * recip_support / c["n_support"], 1),
            "still_fixed": still_fixed,
            "samples_failing": ";".join(failing),
        })

    with open(args.out, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=list(rows[0].keys()) if rows else [
            "chrom", "start", "end", "cnv_type", "claimed_support",
            "reciprocal_support", "reciprocal_pct_of_claimed", "still_fixed",
            "samples_failing"])
        w.writeheader()
        for r in rows:
            w.writerow(r)

    print(f"\n{n_downgraded} of {len(fixed)} 'fixed' CNVRs FAIL strict "
          f">= {args.min_reciprocal_pct}% reciprocal overlap in at least "
          f"one sample and should be re-labeled.")
    print(f"Full report: {args.out}")
    print()
    print("If BMPR1B or the Chr.1 NK-cytotoxicity cluster CNVRs are among")
    print("the downgraded rows, flag that explicitly -- those are the two")
    print("headline findings in the Abstract/Highlights and would need")
    print("re-wording (e.g. 'segregating, 19/21' rather than 'fixed').")
    print("If they are NOT among the downgraded rows, that is itself a")
    print("useful robustness statement to add to the Discussion.")


if __name__ == "__main__":
    main()
