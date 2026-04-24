#!/bin/bash

echo -e "chrom\tpos\tref\talt\tsv_type\tsv_length\tnumber_svs_merged\t$(bcftools query -l ../../.processed_data/structural_variants/141_SVs_merged.vcf | \
	sed -E 's/^[0-9]+_//' | \
	paste -sd '\t' -)" > ../../processed_data/structural_varaints/Jasmine_merged_SVs.tsv

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%SVTYPE\t%SVLEN\t%SUPP[\t%GT]\n' ../../processed_data/structural_varaints/141_SVs_merged.vcf >> ../../processed_data/structural_varaints/Jasmine_merged_SVs.tsv
