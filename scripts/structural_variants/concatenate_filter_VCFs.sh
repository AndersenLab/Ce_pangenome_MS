#!/bin/bash

for file in *.vcf; do bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t%INFO/SVTYPE\t%INFO/SVLEN' $file | awk -v strain=${file%%.*} -v OFS='\t' '$7 >= 50 || $7 <= -50 {print $1,$2,$3,$4,$5,$6,$7,strain}'; done | grep -w "PASS" | grep -v -w "SNV"
