#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(regioneR))

message("Loading CNVRs and Genomic Features...")
cnv_df <- read.table("./02_cnvrs/final_cnvrs.bed", sep="\t", header=FALSE, stringsAsFactors=FALSE)
gff_df <- read.table("./03_annotation/genes_only.gff", sep="\t", comment.char="#", stringsAsFactors=FALSE, quote="")

# 1. Dynamically find the exact chromosomes shared between your CNVRs and the Gene annotations
shared_chrs <- intersect(unique(cnv_df$V1), unique(gff_df$V1))
message(sprintf("Found %d matching chromosomes/scaffolds between datasets.", length(shared_chrs)))

# Filter out mismatched contigs
cnv_df <- cnv_df[cnv_df$V1 %in% shared_chrs, ]
gff_df <- gff_df[gff_df$V1 %in% shared_chrs, ]

# Convert to GenomicRanges
cnvrs <- toGRanges(cnv_df[, 1:4])
features <- toGRanges(data.frame(chr=gff_df$V1, start=as.numeric(gff_df$V4), end=as.numeric(gff_df$V5)))

# 2. Build the reference genome dynamically from only the shared chromosomes
goat_genome <- data.frame(
  chr = shared_chrs,
  start = 1,
  end = 150000000 
)

message("Running Permutation Test (100 iterations for speed)...")
# 3. Disable force.parallel to prevent silent core crashes on Ubuntu
pt <- overlapPermTest(A = cnvrs, B = features, ntimes = 100, genome = goat_genome, force.parallel = FALSE)

print(summary(pt))

pdf("./03_annotation/regioneR_permutation_plot.pdf")
plot(pt)
dev.off()
message("Permutation test complete. Plot saved to ./03_annotation/regioneR_permutation_plot.pdf")
