# bbgoatchip-cnv-analysis

Analysis code for:

> Shil, S. et al. Genome-wide copy number variation landscape of the Black
> Bengal Goat reveals immune-pathway enrichment and QTL co-location.
> *[Journal, in preparation]*.

This repository contains the scripts used to call, filter, and characterize
copy number variants (CNVs) from whole-genome resequencing data of 21
quality-controlled Black Bengal Goat (*Capra hircus*) individuals, and to
run the downstream permutation and functional-enrichment analyses reported
in the manuscript. It does not contain sequencing data or analysis outputs;
those are deposited separately on Zenodo (see below).

## Data availability

The CNVR catalogue, per-sample calls, QC filtering tables, bin-size
sensitivity results, and enrichment/permutation outputs produced by these
scripts are archived at:

Zenodo: https://doi.org/[ZENODO_DOI]

Raw whole-genome resequencing data (BAM/FASTQ) are not yet publicly
available; see the manuscript's Data Availability Statement.

## Repository structure

```
cnv_calling/
  01_cnvpytor_pipeline.sh          rd -> his -> partition -> call, per sample, 100 bp bins
  02_bin_size_sensitivity.sh       re-bins existing rd data at 500/1000 bp for a 6-sample subsample
  03_filter_and_report.py          post hoc per-call quality filter (e-val1, q0, pN, dG) + funnel report

cnvr_construction/
  build_population_cnvrs.py        merges per-sample calls into population CNVRs, singleton removal
  reciprocal_overlap_check.py      per-sample >=50% both-direction overlap confirmation

annotation/
  cnvr_gene_overlap.py             CNVR-gene intersection (BEDTools wrapper)
  cnvr_qtl_overlap.py              CNVR-QTL intersection against AnimalQTLdb

enrichment/
  run_gprofiler.py                 g:Profiler enrichment queries
  run_shinygo.R                    ShinyGO KEGG enrichment

permutation/
  gene_overlap_permutation.py      1,000-iteration null-model permutation, gene overlap
  qtl_overlap_permutation.py       1,000-iteration null-model permutation, QTL overlap

manta_validation/
  run_manta.sh                     structural-variant validation at candidate loci

env/
  environment.yml                  conda environment (cnvpytor, pysam, bedtools, etc.)
```

## Requirements

- CNVpytor
- Python 3.9+ (pysam, pandas, numpy, scipy)
- BEDTools
- R (for ShinyGO enrichment step)
- See `env/environment.yml` for the full pinned environment.

## Reference genome

*Capra hircus* ARS1 reference genome assembly (GCF_001704415.2, ARS1.2).

## Usage

Scripts are numbered/grouped by pipeline stage and intended to be run in
the order listed above. Each script accepts `--help` for usage details.
Paths to input BAM files and output directories are set via a config file
(see individual script headers) and are not included in this repository.

## Citation

If you use this code, please cite the manuscript above.

## License

MIT License — see `LICENSE`.

## Contact

Dr. Sharadindu Shil
Genomics and Molecular Biology Laboratory (GMBL)
Paschimbanga Go Sampad Bikash Sanstha (PBGSBS)
Animal Resources Development Department, Government of West Bengal
dr.sharadindu@gmail.com
ORCID: 0000-0002-6167-2244
