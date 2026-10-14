#!/bin/bash

#SBATCH -J RNA_cov
#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 12:00:00
#SBATCH -c 2
#SBATCH --mem=16G

set -euo pipefail

#GFF_LIST="/vast/eande106/projects/Lance/THESIS_WORK/misc/list_of_GFFs.txt"

#GFF=$(sed -n ${SLURM_ARRAY_TASK_ID}p $GFF_LIST)

GFF="/vast/eande106/data/c_elegans/genomes/PRJNA13758/WS283/csq/c_elegans.PRJNA13758.WS283.csq.gff3"
BAM="/vast/eande106/projects/Lance/THESIS_WORK/manuscript_repos/pangenome_Ce_MS/processed_data/gwas/rna_aln/N2/N2_Aligned.sortedByCoord.out.bam"
strain="N2"

#base=$(basename $GFF)
#strain=${base%%.*}

#RNA_path="/vast/eande106/projects/Lance/THESIS_WORK/gene_annotation/processed_data/rna_WS_genome_aln/STAR_output"
OUT="/vast/eande106/projects/Lance/THESIS_WORK/manuscript_repos/pangenome_Ce_MS/processed_data/gwas/rna_aln/coverage"

#BAM="$RNA_path/$strain/${strain}_Aligned.sortedByCoord.out.bam"

mkdir -p $OUT

if [[ -f $BAM ]]; then
    echo "Calculating coverage for $strain"

    samtools bedcov \
        <(awk -F'\t' 'BEGIN {OFS="\t"} 
    	    $3 == "gene" {
                gene = $9
                sub(/.*ID=/, "", gene)
                sub(/;.*/, "", gene)
                print $1, $4 - 1, $5, gene
            }' $GFF) $BAM |
    awk 'BEGIN {OFS="\t"} {
        region_len = $3 - $2
        mean = (region_len > 0 ? $NF / region_len : 0)
        print $0, mean}' > $OUT/${strain}_rna_gene_coverage.txt
else
    echo "No RNA-seq alignment for $strain, skipping..."
    exit 1
fi
