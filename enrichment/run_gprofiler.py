from gprofiler import GProfiler
import pandas as pd
from datetime import datetime

gp = GProfiler(return_dataframe=True)

with open("fixed_genes_for_enrichment.txt") as f:
    genes = [l.strip() for l in f if l.strip()]

results = gp.profile(
    organism="chircus",
    query=genes,
    sources=["KEGG", "GO:BP", "GO:MF", "GO:CC"],
    user_threshold=0.05,
    significance_threshold_method="fdr",
)

ts = datetime.now().strftime("%Y%m%d_%H%M%S")
results.to_csv(f"gprofiler_fixed_254_{ts}.csv", index=False)
print(results.head(10))
