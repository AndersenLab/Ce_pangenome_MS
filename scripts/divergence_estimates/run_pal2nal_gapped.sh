#!/bin/bash

#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 4:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=pal2nal_array

#source activate pal2nal

#MSA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../processed_data/genome_resources/annotation/single_copy_msa_orthogroups/msa_files.txt)

#BASE=$(basename $MSA .aligned.fa)

#CDS="../../processed_data/divergence_estimates/CDS_fasta/orthogroups/${BASE}.cds.fa"
#OUT="../../processed_data/divergence_estimates/CDS_fasta/codon_aware_msa_wgaps/${BASE}.aligned.codons.gapped.fa"

#if [[ ! -f $CDS ]]; then
 # echo "Missing CDS file: $CDS" >&2
 # exit 1
#fi

#pal2nal.pl $MSA $CDS -output fasta > $OUT

pal2nal.pl ../../processed_data/divergence_estimates/single_copy_msa_orthogroups/OG0003616.fa ../../processed_data/divergence_estimates/CDS_fasta/orthogroups/OG0003616.cds.fa -output fasta > ../../processed_data/divergence_estimates/codon_alignments/OG0003616.aligned.codons.gapped.fa
