#!/bin/bash

module load samtools/1.8

cd ../data

# To obtain bam files sorted and filtered for alignments with a mapping quality (MAPQ) >= 20

parallel -j "4" samtools view {}_sort.bam -b -q 20 ">" {}_filtered.bam ::: $(ls -1 *_sort.bam | sed 's/_sort.bam//')
