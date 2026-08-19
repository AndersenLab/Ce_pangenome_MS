#!/bin/bash

#SBATCH -A eande106
#SBATCH -p parallel
#SBATCH -t 12:00:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --job-name=pal2nal_array

source activate pal2nal

#MSA=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ../../processed_data/divergence_estimates/single_copy_msa_orthogroups/msa_files.txt)

#BASE=$(basename $MSA .aligned.fa)

#CDS="../../processed_data/divergence_estimates/CDS_fasta/orthogroups/${BASE}.cds.fa"
#OUT="../../processed_data/divergence_estimates/CDS_fasta/codon_aware_msa_wgaps/${BASE}.aligned.codons.gapped.fa"

#if [[ ! -f $CDS ]]; then
#  echo "Missing CDS file: $CDS" >&2
#  exit 1
#fi

#pal2nal.pl $MSA $CDS -output fasta > $OUT

for file in ../../processed_data/divergence_estimates/single_copy_msa_orthogroups/renamed_MSA/*MSA.fa; do
	BASE=$(basename $file _seqID_fixed.MSA.fa)
	echo "=========== PROCESSING $BASE ===========" 

	CDS="../../processed_data/divergence_estimates/CDS_fasta/orthogroups/reordered_cds_seqIDs/${BASE}.cds_reordered.fa"
	OUT="../../processed_data/divergence_estimates/codon_alignments/${BASE}.aligned.codons.gapped.fa"

	pal2nal.pl $file $CDS -output fasta > $OUT
done
