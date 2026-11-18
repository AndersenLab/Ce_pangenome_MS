#!/bin/bash

set -euo pipefail

GFF=$(sed -n ${SLURM_ARRAY_TASK_ID}p $GFF_LIST)

base=$(basename $GFF)
strain=${base%%.*}

RNA_path="../../processed_data/gwas/rna_aln"
OUT="../../processed_data/gwas/rna_aln/coverage"

BAM="$RNA_path/$strain/${strain}_Aligned.sortedByCoord.out.bam"

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
