#!/bin/bash

input="list_of_all_PAV_vcfs.tsv"

jasmine --output_genotypes threads=48 outdir=$TMP file_list=$input out_file=../../processed_data/structural_variants/141_SVs_merged.vcf
