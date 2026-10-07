#!/bin/bash

source activate star

base=$(basename $file)
strain=${base%%.*}

if [[ ! -f $RNA_dir/$strain/${strain}.merged_r1.fq.gz ]]; then
    echo "RNA files not found for $strain, skipping"
    exit 0
fi

echo "Aligning RNA to genome for strain $strain ..."

mkdir -p $OUT/STAR_genome_index/$strain
mkdir -p $OUT/STAR_output/$strain

# Genome indexing
STAR \
--runThreadN 24 \
--runMode genomeGenerate \
--limitGenomeGenerateRAM 600000000000 \
--genomeDir $OUT/STAR_genome_index/$strain \
--genomeFastaFiles $file \
--genomeSAindexNbases 12

# Alignment
STAR \
--runThreadN 24 \
--genomeDir $OUT/STAR_genome_index/$strain \
--outSAMtype BAM SortedByCoordinate \
--twopassMode Basic \
--readFilesCommand zcat \
--alignIntronMax 10000 \
--limitBAMsortRAM 5000000000 \
--outSAMstrandField intronMotif \
--outFileNamePrefix $OUT/STAR_output/$strain/${strain}_ \
--readFilesIn $RNA_dir/$strain/${strain}.merged_r1.fq.gz $RNA_dir/$strain/${strain}.merged_r2.fq.gz

echo "Done with strain $strain"
