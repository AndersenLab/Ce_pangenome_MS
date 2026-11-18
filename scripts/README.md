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
	- This script creates Supplementary Figure 26.
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
	- This script looks at quantifying genes that are private/common among wild strains but missing in the N2 reference strain, and vice versa. 
	- This script creates Supplementary Figures 34 and 35.
- `genome_stats.R`
	- This script is used to visualize genome quality statistics.
	- This script creates Supplementary Figure 2
- `haplotypePlotter.R`
	- This script uses genome-genome alignments and orhtology to visualize haplotypes.
	- This script creates Supplementary Figure 15.
- `InterProScan.sh`
	- This script runs InterProScan (IPR) to annotates functional domains of proteomes.
- `nucmer_WS_aln.sh`
	- This script performs genome-genome alignments between a given wild strain and the N2 reference genome using NUCmer.

## gwas/

This directory contains scripts used to analyze data pertaining to the abamectin resistance genome-wide association study.

### Contains

- `abamectin_manhattan_plot.R`
	- This script is used to visualize SNVs that are associated with differences in response to abamectin.
	- This script creates Supplementary Figure 27.
- `haplotype_plotter_abamectin_ivermectin_20260922.R`
	- This script is used to visualize haplotype diversity of wild strains surrounding the peak marker associated with abamectin resistance.
	- This script creates Figure 6.
	- This script creates Supplementary Figures 28 and 29.
- `rna_coverage_stat.sh`
	- This script calculates RNA alignment coverage of genes for every wild strain and the N2 reference strain.
- `rna_expression.R`
	- This script calculates stats for RNA fold-coverage of genes in the QTL associated with abamectin resistance.
	- This script creates Supplementary Figure 31.
- `rna_ws_genome_aln.sh`
	- This script aligns RNA to wild strain and N2 genoems using STAR.
- `snvs_in_interval.R`
	- This script gathers statistics and data on SNVs found in the QTl associated with abamectin resistance.
	- This script creates Supplementary Figure 30.
	- This script creates Supplementary Table 15.

## hyperdivergent_regions/

This directory performs HDR liftover to wild strain genomic coordinates using genome-genome alignments and perform statistics on the types of genes found in wild strains HDRs and their conservation with reference strain HDR genes.

### Contains
- `cnv_pav_HDR_enrichment.R`
	- This looks at wild strain HDR gene content copy-number conservation with N2 reference strain HDR genes.
	- This script creates Figure 5.
	- This script creates Supplementary Figures 22, 23, 24, and 25.
- `hdr_liftover.R`
	- This script peforms HDR coordinate lift-over using genome-genome alignments with the N2 reference genome.
	- This script creates Supplementary Figures 16, 17, and 18.
	- This script creates Supplementary Tables 6, 7, and 8.
- `wild_strain_HDR_enrichment.R`
	- This script looks at HDR gene content in wild strains and performs functional enrichment analyses to understand if certain gene classes are over-represented in HDRs.
	- This script creates Figure 4.
	- This script creates Supplementary Figures 19, 20, and 21.
	- This script creates Supplementary Table 10.

## orthology/

This directory

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


