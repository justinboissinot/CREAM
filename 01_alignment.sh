#!/bin/bash

module load sabre/1.000
module load cutadapt/3.2
module load bwa/0.7.17
module load samtools/1.8


cd ../data


# 1. Demultiplexing of the paired-end reads with Sabre

parallel -j "4" sabre pe -f {}_R1.fq.gz -r {}_R2.fq.gz -b {}_PE -u {}_R1.unknown.barcodes -w {}_R2.unknown.barcodes ::: $(ls -1 *_R1.fq.gz | sed 's/_R1.fq.gz//')


# 2. Trimming of the reads with cutadapt

adapfor=AGATCGGAAGA
adaprev=NNNNNNNAGATCGGAAGA # Since sabre only removes barcodes on one end, we added 7 N bases to fully remove even the longest barcodes 
readlen=50

parallel -j "4" cutadapt -O 10 -a ${adapfor} -A ${adaprev} -m ${readlen} -o {}_R1.fastq -p {}_R2.fastq {}_R1.fq {}_R2.fq ::: $(ls -1 *_R1.fq | sed 's/_R1.fq//')


# 3. Alignment of trimmed reads to the reference genome with BWA-MEM

refgen=GCF_900626175.2_cs10_genomic.fna

parallel -j "4" bwa mem -t "2" ../refgenome/"${refgen}" {}_R1.fastq {}_R2.fastq ">" {}.sam ::: $(ls -1 *_R1.fastq | sed 's/_R1.fastq//')

parallel -j "4" samtools view -b -S -h {}.sam ">" {}_temp.bam ::: $(ls -1 *.sam | sed 's/.sam//')

parallel -j "4" samtools sort {}_temp.bam -o {}_sort.bam ::: $(ls -1 *_temp.bam | sed 's/_temp.bam//')

rm *temp.bam
rm *.sam
