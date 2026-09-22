#!/usr/bin/env python3
"""
06_run_corrected_gprofiler.py

Runs the corrected (246-CNVR, reciprocal-overlap-confirmed) gene list from
05_rebuild_fixed_enrichment.py through g:Profiler, using IDENTICAL parameters
to the original /home/gmbluser/run_gprofiler.py (organism=chircus,
sources=KEGG/GO:BP/GO:MF/GO:CC, FDR<0.05) so the result is directly
comparable to the original 254-gene, p=1.5e-21 NK-cytotoxicity finding.

USAGE:
  python3 06_run_corrected_gprofiler.py corrected_fixed_genes_246cnvr.txt
"""
import sys
from datetime import datetime
from gprofiler import GProfiler
import pandas as pd

if len(sys.argv) != 2:
    sys.exit("Usage: python3 06_run_corrected_gprofiler.py <gene_list.txt>")

gene_list_path = sys.argv[1]
gp = GProfiler(return_dataframe=True)

with open(gene_list_path) as f:
    genes = [l.strip() for l in f if l.strip()]

print(f"Running g:Profiler on {len(genes)} genes from {gene_list_path} ...")

results = gp.profile(
    organism="chircus",
    query=genes,
    sources=["KEGG", "GO:BP", "GO:MF", "GO:CC"],
    user_threshold=0.05,
    significance_threshold_method="fdr",
)

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
out_path = f"corrected_fixed_enrichment_{ts}.csv"
results.to_csv(out_path, index=False)
print(f"\nFull results written to: {out_path}  ({len(results)} significant terms)")

pd.set_option("display.max_colwidth", 60)

nk = results[results["native"] == "KEGG:04650"]
print("\n=== NK-cell mediated cytotoxicity (KEGG:04650) — the original")
print("=== headline result. Empty output below means it did NOT survive")
print("=== in the corrected, strictly-confirmed-fixed gene set:")
print(nk[["source", "native", "name", "p_value", "intersection_size", "term_size"]]
      if not nk.empty else "  (no rows -- NK-cytotoxicity is NOT significant in the corrected set)")

print("\n=== Top 15 terms overall in the corrected result, by p-value:")
cols = [c for c in ["source", "native", "name", "p_value", "intersection_size", "term_size"] if c in results.columns]
print(results.sort_values("p_value")[cols].head(15).to_string(index=False))
