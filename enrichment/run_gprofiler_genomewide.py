from gprofiler import GProfiler
import pandas as pd
from datetime import datetime

gp = GProfiler(return_dataframe=True)

with open("unique_cnv_genes_clean.txt") as f:
    genes = [l.strip() for l in f if l.strip()]

print(f"Loaded {len(genes)} genes")  # sanity check — should print 2856

results = gp.profile(
    organism="chircus",
    query=genes,
    sources=["KEGG"],
    user_threshold=0.05,
    significance_threshold_method="fdr",
)

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
results.to_csv(f"gprofiler_genomewide_2856_{ts}.csv", index=False)
print(f"Total significant KEGG pathways: {len(results)}")
print(results[["native","name","p_value","term_size","query_size","intersection_size"]].head(10).to_string())
