#!/bin/bash

# align with nucmer (will spit out a .delta file)
strain=$(basename $file .fa)
nucmer --maxgap=500 --prefix=${strain}.contigs --coords ../../processed_data/genome_resources/genome_data/c_elegans.PRJNA13758.WS283.genome.fa ${file}

# get coordinate file - filter to alignments that are >1kb
show-coords -r -l -T ${strain}.contigs.delta | awk '$5 > 1000' > ../../processed_data/genome_resources/genome_data/${strain}.transformed.tsv
