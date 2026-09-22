import re
from collections import defaultdict

# ---- 1. Parse gene annotations (marker_id -> gene) from earlier step ----
marker_gene = {}
with open("FINAL_validated_genes.tsv") as f:
    next(f)  # header
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 2:
            continue
        marker_gene[cols[0].strip()] = cols[1].strip()

# ---- 2. Parse nearest-QTL results, keeping ALL ties at the minimum distance ----
# Columns: chrom,start,end,marker_id,svtype (5) | qtl_chrom,qtl_start,qtl_end,qtl_name (4) | distance (1)
marker_best_dist = {}
marker_qtls = defaultdict(set)
marker_svtype = {}

with open("elite_nearest_qtl.txt") as f:
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 10:
            continue
        marker_id = cols[3]
        svtype = cols[4]
        qtl_name = cols[9-1]  # col 9 (0-indexed 8) = QTL name field
        dist = int(cols[9])   # col 10 (0-indexed 9) = distance

        marker_svtype[marker_id] = svtype

        if marker_id not in marker_best_dist or dist < marker_best_dist[marker_id]:
            marker_best_dist[marker_id] = dist
            marker_qtls[marker_id] = {qtl_name}
        elif dist == marker_best_dist[marker_id]:
            marker_qtls[marker_id].add(qtl_name)

# ---- 3. Build final table ----
THRESHOLD = 250000  # bp -- adjust if needed

rows = []
for marker_id, dist in marker_best_dist.items():
    svtype = marker_svtype[marker_id]
    gene = marker_gene.get(marker_id, "NA")
    qtls = "; ".join(sorted(marker_qtls[marker_id]))
    within_threshold = "YES" if dist <= THRESHOLD else "no"

    # Tier logic: top tier = DEL-type (statistically validated), has a gene, AND near a QTL
    if svtype == "DEL" and gene != "NA" and dist <= THRESHOLD:
        tier = "Tier 1: DEL + genic + QTL-proximal"
    elif gene != "NA" and dist <= THRESHOLD:
        tier = "Tier 2: genic + QTL-proximal"
    elif dist <= THRESHOLD:
        tier = "Tier 3: QTL-proximal only"
    elif gene != "NA":
        tier = "Tier 4: genic, no nearby QTL"
    else:
        tier = "Tier 5: neither"

    rows.append((marker_id, svtype, gene, qtls, dist, within_threshold, tier))

# ---- 4. Sort by tier, then distance ----
tier_order = {t: i for i, t in enumerate([
    "Tier 1: DEL + genic + QTL-proximal",
    "Tier 2: genic + QTL-proximal",
    "Tier 3: QTL-proximal only",
    "Tier 4: genic, no nearby QTL",
    "Tier 5: neither",
])}
rows.sort(key=lambda r: (tier_order[r[6]], r[4]))

with open("QTL_GENE_TARGETING_FINAL.tsv", "w") as out:
    out.write("marker_id\tsvtype\tgene\tnearest_QTL(s)\tdistance_bp\twithin_250kb\ttier\n")
    for r in rows:
        out.write("\t".join(str(x) for x in r) + "\n")

# ---- 5. Summary ----
from collections import Counter
tier_counts = Counter(r[6] for r in rows)
print(f"Total markers: {len(rows)}")
for t in tier_order:
    print(f"  {t}: {tier_counts.get(t, 0)}")

print("\nTop 15 Tier 1 candidates:")
tier1 = [r for r in rows if r[6] == "Tier 1: DEL + genic + QTL-proximal"]
for r in tier1[:15]:
    print(f"  {r[0]}  gene={r[2]}  dist={r[4]}bp  QTL={r[3][:80]}")

print(f"\nSaved: QTL_GENE_TARGETING_FINAL.tsv")
