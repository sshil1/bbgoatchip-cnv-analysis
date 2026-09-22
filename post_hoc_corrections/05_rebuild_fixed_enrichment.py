#!/usr/bin/env python3
"""
05_rebuild_fixed_enrichment.py

PURPOSE
  Re-derives the manuscript's headline "fixed CNVR" gene-set enrichment
  (originally: 254 genes, NK-cytotoxicity top hit, p=1.5e-21) using only the
  246 CNVRs that survive STRICT reciprocal overlap (>=50%, script 04), instead
  of all 1,272 CNVRs that were merely proximity-merged with 21/21 sample
  support. This directly answers whether the NK-cytotoxicity finding is an
  artifact of the loose "fixed" definition or a robust result.

  Context (established 2026-09-18): script 04 showed only 246/1,272 CNVRs
  labeled "fixed" in population_cnvrs_filtered_clean.bed actually have all 21
  samples' individual calls reciprocally overlapping >=50%. The Chr.1
  NK-cytotoxicity cluster region specifically FAILED strict reciprocal
  overlap in all 5 of its CNVRs checked -- this script determines whether the
  broader 254-gene enrichment result survives on the corrected, stricter set.

  The original intermediate script that built fixed_genes_for_enrichment.txt
  from population_cnvrs_filtered_clean.bed could not be located on the
  server (not a reproducibility red flag on its own -- it was one script
  among many from a fast-moving analysis -- but it means this rebuild uses a
  fresh, clean intersection rather than literally re-running the original
  code). It uses the SAME gene annotation file referenced elsewhere in this
  project (goat_ARS1.2_genes_only.gff.gz, used by rank_sv_markers.py) and the
  SAME g:Profiler call parameters as run_gprofiler.py (organism=chircus,
  sources=KEGG/GO:BP/GO:MF/GO:CC, FDR<0.05), so the comparison to the
  original 254-gene result is apples-to-apples on methodology, differing
  only in which CNVR set feeds it.

INPUTS
  --regenotype-report   output of 04_fixed_cnvr_regenotype_check.py, run
                         against the UNFILTERED calls (the correct one to use
                         -- see 02026-09-18 conversation notes). Real path:
                         ~/cnv_repro_scripts/fixed_cnvr_regenotype_report_unfiltered.tsv
  --gff                 gene annotation, gzipped OK. Real path:
                         /home/gmbluser/goat_ARS1.2_genes_only.gff.gz
  --original-gene-list  fixed_genes_for_enrichment.txt, for side-by-side
                         reporting of which genes were lost/kept (optional
                         but recommended). Real path:
                         /home/gmbluser/fixed_genes_for_enrichment.txt

OUTPUT
  corrected_fixed_genes_246cnvr.txt  -- new gene list, one symbol per line
  corrected_fixed_enrichment_<timestamp>.csv -- g:Profiler result
  Console summary: gene-count comparison + whether NK-cytotoxicity
  (KEGG:04650) survives in the corrected result, with its new p-value.

USAGE (run on gmbluser@10.0.0.2, conda env with `gprofiler-official` and
  internet access -- same env run_gprofiler.py used):
  python3 05_rebuild_fixed_enrichment.py \
      --regenotype-report ~/cnv_repro_scripts/fixed_cnvr_regenotype_report_unfiltered.tsv \
      --gff /home/gmbluser/goat_ARS1.2_genes_only.gff.gz \
      --original-gene-list /home/gmbluser/fixed_genes_for_enrichment.txt
"""
import argparse
import csv
import gzip
import sys
from collections import defaultdict
from datetime import datetime


def open_maybe_gz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_genes(gff_path):
    """Same parsing logic as rank_sv_markers.py's load_genes(), for
    consistency with the rest of this project's pipeline."""
    genes = defaultdict(list)
    n = 0
    with open_maybe_gz(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            chrom, start, end, attrs = f[0], int(f[3]), int(f[4]), f[8]
            name = None
            for kv in attrs.split(";"):
                kv = kv.strip()
                if kv.startswith("Name=") or kv.startswith("gene="):
                    name = kv.split("=", 1)[1]
                    break
                if kv.startswith("ID=gene-"):
                    name = kv.split("ID=gene-", 1)[1]
            genes[chrom].append((start, end, name or "unnamed"))
            n += 1
    print(f"Loaded {n} gene features across {len(genes)} chromosomes from {gff_path}",
          file=sys.stderr)
    return genes


def load_confirmed_fixed_cnvrs(report_path):
    """Reads script 04's report (comma-delimited), returns list of
    (chrom, start, end) for rows where still_fixed == 'True'."""
    cnvrs = []
    with open(report_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("still_fixed") == "True":
                cnvrs.append((row["chrom"], int(row["start"]), int(row["end"])))
    return cnvrs


def overlapping_genes(chrom, start, end, gene_index):
    hits = set()
    for gstart, gend, gname in gene_index.get(chrom, []):
        if start <= gend and end >= gstart:
            hits.add(gname)
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--regenotype-report", required=True)
    ap.add_argument("--gff", required=True)
    ap.add_argument("--original-gene-list", default=None)
    ap.add_argument("--out-gene-list", default="corrected_fixed_genes_246cnvr.txt")
    args = ap.parse_args()

    confirmed = load_confirmed_fixed_cnvrs(args.regenotype_report)
    print(f"Loaded {len(confirmed)} strictly-confirmed-fixed CNVRs from "
          f"{args.regenotype_report}", file=sys.stderr)
    if len(confirmed) == 0:
        sys.exit("No rows with still_fixed=True found -- check the report "
                  "file / column name.")

    gene_index = load_genes(args.gff)

    corrected_genes = set()
    for chrom, start, end in confirmed:
        corrected_genes |= overlapping_genes(chrom, start, end, gene_index)
    corrected_genes.discard("unnamed")

    with open(args.out_gene_list, "w") as out:
        for g in sorted(corrected_genes):
            out.write(g + "\n")

    print(f"\nCorrected (246-CNVR, reciprocal-overlap-confirmed) fixed gene "
          f"set: {len(corrected_genes)} genes -> {args.out_gene_list}")

    if args.original_gene_list:
        with open(args.original_gene_list) as fh:
            original_genes = set(l.strip() for l in fh if l.strip())
        kept = corrected_genes & original_genes
        lost = original_genes - corrected_genes
        gained = corrected_genes - original_genes
        print(f"Original (1,272-CNVR) fixed gene set: {len(original_genes)} genes")
        print(f"  Kept in corrected set:  {len(kept)}")
        print(f"  Lost (were in original, not in corrected): {len(lost)}")
        print(f"  Gained (new, not in original): {len(gained)}")
        for marker in ["CD48", "CD244", "C3"]:
            status = "KEPT" if marker in corrected_genes else (
                "LOST" if marker in original_genes else "not in original either")
            print(f"  Marker gene {marker}: {status}")

    print("\nNext: run this gene list through g:Profiler with the SAME "
          "parameters as run_gprofiler.py:")
    print(f"""
from gprofiler import GProfiler
import pandas as pd
from datetime import datetime

gp = GProfiler(return_dataframe=True)
with open("{args.out_gene_list}") as f:
    genes = [l.strip() for l in f if l.strip()]

results = gp.profile(
    organism="chircus",
    query=genes,
    sources=["KEGG", "GO:BP", "GO:MF", "GO:CC"],
    user_threshold=0.05,
    significance_threshold_method="fdr",
)
ts = datetime.now().strftime("%Y%m%d_%H%M%S")
results.to_csv(f"corrected_fixed_enrichment_{{ts}}.csv", index=False)
print(results[results['native'] == 'KEGG:04650'])  # NK-cytotoxicity, if present
print(results.sort_values('p_value').head(15))
""")


if __name__ == "__main__":
    main()
