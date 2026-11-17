## Scripts
This directory contains the scripts that were used to perform analyses and generate figures and tables for the manuscript.

## Directory Structure

```
script/
├── divergence_estimates
├── genomes
├── gwas
├── hyperdivergent_regions
├── orthology
└── structural_variants
```

## divergence_estimates/

This directory contains scripts to calculate Ka and Ks for single-copy orthologs.

### Contains

- `get_cds_fasta.sh`
	- This script creates codon FASTA files from GFF files.
- `ka_ks_calculations.R`
	- This script calcuales Ka and Ks for single-copy orthologs.
	- This script creates Supplementary Figure 28
- `protein_msa.sh`
	- Multiple sequence alignments of single-copy orthologs.
- `run_pal2nal_gapped.sh`
	- This script runs PAL2NAL on single-copy orthologs to create codon MSAs.

## genomes/

This directory contains scripts to visualize/calculate genome stats, haplotypes, genome-genome alignments, and perform InterProScan annotation.

### Contains
- `elegans_map_relatedness_asmStats.R`
	- Concatenates world map indicating where wild strains have been collected, relatedness matrix, and genome assembly stats to create Figure 1
