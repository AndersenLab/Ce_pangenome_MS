#!/bin/bash

#SBATCH -J cds_fasta
#SBATCH -A eande106          # Allocation name
#SBATCH -p parallel          # Partition/Queue name
#SBATCH -t 12:00:00          # Job walltime/duration (hh:mm:ss)
#SBATCH -N 1                 # Number of nodes
#SBATCH -n 12                 # Number of cores

# Loop through GFFs and produce spliced CDS FASTAs
for GFF in ../../processed_data/genome_resources/annotation/gffs/*.gff3; do
    SAMPLE=$(basename $GFF .longest.gff3)
    GENOME="../../processed_data/genome_resources/genome_data/genomes/${SAMPLE}.genome.fa"
    OUT="../../processed_data/divergence_estimates/CDS_fasta/${SAMPLE}.spliced_cds.fa"

    if [[ -f $GENOME ]]; then
        echo "Processing $SAMPLE"
        gffread -x $OUT -g $GENOME $GFF
    else
        echo "Genome not found for $SAMPLE: $GENOME" >&2
    fi
done
