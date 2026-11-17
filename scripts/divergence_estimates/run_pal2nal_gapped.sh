#!/bin/bash

source activate pal2nal

for file in ../../processed_data/divergence_estimates/single_copy_msa_orthogroups/renamed_MSA/*MSA.fa; do
	BASE=$(basename $file _seqID_fixed.MSA.fa)
	echo "=========== PROCESSING $BASE ===========" 

	CDS="../../processed_data/divergence_estimates/CDS_fasta/orthogroups/reordered_cds_seqIDs/${BASE}.cds_reordered.fa"
	OUT="../../processed_data/divergence_estimates/codon_alignments/${BASE}.aligned.codons.gapped.fa"

	pal2nal.pl $file $CDS -output fasta > $OUT
done
