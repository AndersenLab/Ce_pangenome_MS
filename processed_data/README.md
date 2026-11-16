# Processed Data

This directory contains the processed datasets used for downstream genomic analyses. Each subdirectory corresponds to a major analysis component. The files are derived from upstream pipelines (see `../scripts/`).

## Directory Structure

processed_data/
- divergence_estimates/
	- codon_alignments/
- genome_resources/  
	- annotation/
		- gffs/
	- genome_data/
		- genomes/
	- isotypes/
	- trees/
- gwas/
	- haemonchus_contortus_orthology/
	- rna_aln/
		- coverage/
	- variant_annotation/
- hdr_liftover/
- orthology/
	- orthofinder/
		- orthofinder_output/
		- proteomes/
			- accessory_78/
			- core_64/
		- proteomes_HC_CE/
			- OrthoFinder/
				- Results_Nov03/
					- Orthogroups/
- structural_variants/

## divergence_estimates/

This directory contains data for synonymous (Ks) and non-synonymous (Ka) divergence estimates for single-copy orthologs. 

### Contains
- ``
	- 

## genome_resources/

This directory contains data pertaining to *C. elegans* genomes assembled in this study, gene annotations, relatedness estimates, and wild strain sampling data.

### Contains
- `annotation/`  
	- `140_wild_strains_IPR_CtypeLectins.tsv`
		- C-type lectins defined by InterProScan (IPR) annotations.
	- `140_wild_strains_IPR_cytochromeP450.tsv`
                - Cytochrome P450s defined by IPR annotations.
        - `140_wild_strains_IPR_fBox.tsv`
                - F-box proteins defined by IPR annotations.
        - `140_wild_strains_IPR_gpcrs.tsv`
                - G protein-coupled receptors (GPCRs) defined by IPR annotations.
        - `140_wild_strains_IPR_nhr.tsv`
                - Nuclear hormone receptors (NHRs) defined by IPR annotations.
        - `140Ws_CGC1_longestIsoGenes_BRAKER.tsv` - Large file. Publicly available on AWS.
                - Gene coordinates for 140 wild strains and CGC1 predicted by BRAKER.
        - `140WSs_CGC1_longestIso_mRNAandGenes_BRAKER.tsv` - Large file. Publicly available on AWS.
                - mRNA and gene coordinates for 140 wild strains and CGC1 predicted by BRAKER.
        - `142_strains_longestIso_mRNAandGenes.tsv` - Large file. Publicly available on AWS.
                - mRNA and gene coordiantes for 140 wild strains, N2, and CGC1 predicted by BRAKER.
        - `c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.gff3`
                - Longest isoform gene annotations for N2 WormBase release WS283.
        - `c_elegans.PRJNA13758.WS283.csq.PCfeaturesOnly.longest.w-mRNA_featuresgff3` - Large file. Publicly available on AWS.
                - Longest isoform gene annotations with associated mRNA coordinates for N2 WormBase release WS283.
        - `chromosome_domain_Celegans.csv`
                - Chromsome tip, arm, and center domains defined in Rockman *et al., PLoS Genetics* 2009. 
	- `gffs/`
		- `*.longest.gff3`
			- The longest-isoform gene models for 140 wild strains, CGC1, and N2. Unfiltered GFF files are available on CaeNDR.
	- `haemonchus_contortus.PRJEB506.WBPS19.annotations.gene_only.longest_isoform.gff3` - Large file. Publicly available on AWS.
		- Protein-coding GFF file for *Haemonchus contortus* WormBase ParaSite release WBPS19. 
	- `N2_geneName_seqName.tsv`
		- Gene name and corresponding sequence names for N2 WormBase GFF file release WS283.
	- `WS276.mRNA_coords_tranName.tsv`
		- Transcript information for N2 WormBase GFF file release WS276.
- `genome_data/`
	- `140WSs_biallelicSNPs.tsv`
		- Biallelic single-nucleotide variants among 140 *C. elegans* wild strains in pangenome strain set in relation to the N2 reference genome. 
	- `141_nucmer_ECA741CGC1.tsv`
		- Genome-genome alignments of 140 wild strains and CGC1 against the N2 reference genome.
	- `genomes/`
		- `*.fa`
			- Contiguous *de novo* genome assemblies for 140 *C. elegans* wild strains. These assemblies are publicly availabe on NCBI at PRJNA1558933 and at CaeNDR.
- `isotypes/
	- `684_Ce_concordance.txt`
		- The number of discordant sites (SNVs) and total number of sites among all isotype reference strain pairs. These data are used to create the relatedness heatmap in Figure 1b. 
	- `Ce_isotype_groups.tsv`
		- Isotype reference strains for all isotype groups among all wild strains of *C. elegans*.
	- `elegans_isotypes_sampling_geo_hawaii_islands.tsv`
		- Island resolution for strains collected from Hawaii.
	- `elegans_isotypes_sampling_geo.tsv`
		- Sampling location of all isotype reference strains for *C. elegans*.
- `trees`
	- `BUSCO_supermatrix.fa.contree`
		- Maximum-likelihood tree created from MSA of single-copy orthologs (BUSCOs) shared among 142 strains of *C. elegans*.
	- `Ce_isotypes_LD_0.9.phy.contree`
		- Maximum-likelihood tree created from LD-pruned SNVs of 684 isotype reference strains. This tree is used to create Supplementary Figure 1.

## gwas/

This directory contains data pertaining to the genome-wide association study of resistance to the anthelmintic drug abamectin.

### Contains
- ``
	- 

## hdr_liftover/

This directory contains wild strain HDR liftover coordinates and related statistics and meta data.

### Contains
- ``
	- 

## orthology/

This directory contains output from OrthoFinder and data related to pangenome gene set analyses. 

### Contains
- ``
	- 

## structural_variants/

This directory contains SV calls and related statistics and meta data.

### Contains
- ``
	- 
