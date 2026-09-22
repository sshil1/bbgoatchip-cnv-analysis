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

Zenodo: https://doi.org/\[ZENODO_DOI\]

GitHub: https://github.com/sshil1/bbgoatchip-cnv-analysis

Raw whole-genome resequencing data (BAM/FASTQ) are not yet publicly
available; see the manuscript's Data Availability Statement.

## Repository structure
EOF

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
