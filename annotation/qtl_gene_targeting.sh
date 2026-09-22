#!/bin/bash
cd /home/gmbluser

echo "===== Locating QTL reference file ====="
QTL_BED=$(find /home/gmbluser -maxdepth 4 -iname "*QTL*.bed*" 2>/dev/null | head -1)
if [ -z "$QTL_BED" ]; then
    echo "Not found by name. Searching more broadly..."
    QTL_BED=$(find /home/gmbluser -maxdepth 5 -iname "*qtl*" 2>/dev/null | grep -iE "\.(bed|gff|txt|tsv)$" | head -1)
fi
echo "Using: $QTL_BED"
head -3 "$QTL_BED"
echo "Columns:"
awk -F'\t' '{print NF; exit}' "$QTL_BED"

echo
echo "===== Sorting inputs ====="
sort -k1,1 -k2,2n "$QTL_BED" > qtl_regions.sorted.bed
# elite_580_v2.bed already has: chrom, start, end, marker_id, svtype (from earlier step)
sort -k1,1 -k2,2n elite_580_v2.bed > elite_580_v2.sorted.bed

echo
echo "===== Intersecting elite markers with QTL regions ====="
bedtools intersect -a elite_580_v2.sorted.bed -b qtl_regions.sorted.bed -wao > elite_vs_qtl.txt
echo "Markers with QTL overlap:"
awk -F'\t' '$NF > 0' elite_vs_qtl.txt | cut -f4 | sort -u | wc -l
