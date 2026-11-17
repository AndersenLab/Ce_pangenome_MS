## Scripts
This directory contains the scripts that were used to perform analyses and generate figures and tables for the manuscript.

## Directory Structure

```
scripts/
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
	- Concatenates world map indicating where wild strains have been collected, relatedness matrix, and genome assembly stats to create Figure 1.
- `gene_set_QC.R`
- `genome_stats.R`
- `haplotypePlotter.R`
- `InterProScan.sh`
- `nucmer_WS_aln.sh`

## gwas/

This directory...

### Contains

- `abamectin_manhattan_plot.R`
- `haplotype_plotter_abamectin_ivermectin_20260922.R`
- `rna_coverage_stat.sh`
- `rna_expression.R`
- `rna_ws_genome_aln.sh`
- `snvs_in_interval.R`

## hyperdivergent_regions/

This directory...

### Contains
- `cnv_pav_HDR_enrichment.R`
- `hdr_liftover.R`
- `wild_strain_HDR_enrichment.R`

## orthology/

This directory...

### Contains
- `64_core_selectionOrthoFinder.R`
- `CE_HC_QTL_orthology.R`
- `gene_family_heatmaps.R`
- `gene_set_visualization.R`
- `orthofinder3_haemonchus_contortus.sh`
- `orthofinder3.sh`
- `pca_hawaiian_analyses.R`

## structural_variants/

This directory...

### Contains
- `concatenate_filter_VCFs.sh`
- `jasmine_SVmerging.sh`
- `jasmine_TSV_creation.sh`
- `pav.sh`
- `structural_variants_stats.R`
- `structural_variants_visualization.R`
- `visualization_SVbyEye.R`


