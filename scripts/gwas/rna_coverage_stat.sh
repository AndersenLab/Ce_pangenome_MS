#!/bin/bash

#SBATCH -J RNA_aln
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 8:00:00
#SBATCH -N 1
#SBATCH -c 8

GFF=$file
base=$(basename $GFF)
strain=
RNA_path="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/rna_WS_genome_aln/STAR_output"

# Get coverage per gene
bedtools coverage -a <(awk '$3 == "gene"' $GFF) -f -b $RNA_path/$strain/${strain}_Aligned.sortedByCoord.out.bam -mean > gene_coverage.txt
