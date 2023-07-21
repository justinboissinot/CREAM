
# CREAM : Comparative Restriction Enzyme Analysis of Methylation

Workflow suggested for the CREAM approach.

Developped by Eric Normandeau and Justin Boissinot.

### Reference
This methylotyping workflow has been published within this paper:
[Ref]

## About CREAM

This workflow is used to efficently call methylotypes, binary scores for DNA methylation at specific loci, from raw reads obtained with the Comparative Restriction Enzyme Analysis of Methylation (CREAM) method. 
Briefly, it compares fragments from sequencing libraries generated with restriction enzymes of different sensitivity to DNA methylation and establishes the location of methylated positions across the genome.

**Attention:**
The current code is still hard coded for our data and is not maintained for flexible usage yet. Despite these limitations, this code can be used as inspiration or better understanding of the CREAM approach.

![Figure CREAM](https://github.com/justinboissinot/CREAM/assets/56981551/b4802a87-d890-4802-8923-a31d3f6b4bc3)
**Figure 1.** Schematic representation of the Comparative Restriction Enzyme Analysis of Methylation (CREAM) approach. a) The same genomic DNA was used as input for the preparation of two libraries. Digestion of the DNA molecules was performed using a set of restriction enzymes sharing a same restriction site but with different sensitivity to DNA methylation (1). Then, universal adapters and sample-specific barcodes were ligated to digested fragments (2). DNA fragments were size selected (3), amplified (4) and sequenced (5). b) The comparative analysis of sequencing data of both libraries leads to four possibilities based on either presence/absence or the length of the DNA fragments. Shared fragments in both libraries with the same length and location in the genome indicate the absence of DNA methylation (possibility A) while differences in the presence or the length of the fragments between the libraries indicate the existence of DNA methylation (possibilities B, C and D).


## Overview of main steps

1. Alignment of demultiplexed and trimmed reads to the reference genome with BWA-MEM
2. Quality filtering on BAM files
3. Quality metrics extraction from inserts (pair of aligned reads)
4. Inserts filtration based on quality metrics and abundance
5. Loci (overlapping inserts) extraction
6. Methylation status call on loci

### Work inspired by Fast-GBS_v2
The first main step, from demultiplexing of the reads to their alignment to the reference genome, was inspired by a pipeline developed within our institute called [Fast-GBS_v2](https://bitbucket.org/jerlar73/fast-gbs_v2/src/master/).
Some scripts from this pipeline were used to initiate the data analysis and facilitate the beginning of the current workflow.

### Dependencies
In addition to the scripts available, the following modules were used:
- sabre/1.000
- cutadapt/3.2
- bwa/0.7.17
- samtools/1.8
- python/3.5

### Data preparation
##### Directory architecture
The following directories should be created for the analysis and results storage by running the `./makeDir.sh` script:
```
barcodes
data
refgenome
```

##### Sequence data preparation
Raw FASTQ files from both sequencing libraries should be added the ``data`` directory. Paired-end sequencing files should be named as:
```
flowcell_lane#_R1.fq.gz
flowcell_lane#_R2.fq.gz
```

In the case of the data from the paper, the files from 2 paired-end sequencing libraries were named as such:
```
CanHpaII_1_R1.fq.gz
CanHpaII_1_R2.fq.gz
CanMspI_1_R1.fq.gz
CanMspI_1_R2.fq.gz
```

##### Barcodes preparation
Barcodes should be added the ``barcodes`` directory, with one file per sequencing library. They should have tab-delimited format with the barcode in the first column and the sample name in the second column:
```
TGCA    AA_875
ACTA    AA_914
CAGA    AA_819
...
```
The script ``makeBarcodeSabre.py`` from the Fast-GBS_v2 pipeline was used to create the required barcodes files used by ``Sabre``. The resulting barcodes files, one per sequencing library, are then stored in the ``data`` directory as follows:
```
CanHpaII_1_PE
CanMspI_1_PE
```

##### Reference genome preparation
The reference genome should be added to the ``refgenome`` directory. Only the main FASTA file is needed.

### 01. Alignment of reads to the reference genome

The next steps are ran using ``parallel`` to optimize computational time consumage.

##### Demultiplexing
Raw reads are first demutiplexed with Sabre and assigned to their corresponding sample based on the barcodes. 
Since paired-end reads were generated in the sequencing libraries, the ``pe`` option was used.
Demultiplexed reads are then stored in their new sample's ``.fq`` file.

##### Trimming
Adaptors are removed from the demultplixed reads with cutadapt. 
Trimmed reads are then stored in their new sample's ``.fastq`` file.

##### Alignment to the reference genome
Trimmed reads are aligned to the reference genome using the BWA-MEM algorithm. Resulting SAM files are then sorted and kept as sorted BAM files (``_sort.bam``), with now one file per sample per sequencing library.

### 02. Quality filtering on BAM files
A small script is used to filter the BAM files for alignments with a mapping quality (MAPQ) of over 20. 
This step was added based on an analysis of the MAPQ scores across the samples.

### 03. Quality metrics extraction from inserts (pair of aligned reads)
The filtered BAM files are used as input for the information extraction script. The main steps for this script are as follows:
- Only reads flagged as mapped to their proper pair are kept.
- Quality metrics (information generated by BWA-MEM) are extracted for all reads.
	- CIGAR string: the length of CIGAR strings (indicating insertions or deletions in the sequence)
	- Distance: the number of differences between the sequence and the reference.
	- Score: the number of mismatches in the alignment (alignment score)
- Pairs of reads are formed and called ranges (equivalent to inserts).
- A list of all ranges (with the corresponding quality metrics from the reads that compose them) is outputted in a new file for each sample.

### 04. Inserts filtration based on quality metrics and abundance
The files with the list of ranges for each sample (1 file per sample) are used as input for the the filtration script. The main steps for this script are as follows:
- Ranges, which correspond to lines in the file, are removed if they were aligned on the remaining scaffolds of the reference genome or the mitochondrial genome.
- Inserts are compiled for all samples in each sequencing library (*Msp*I and *Hpa*II). Inserts are here defined as ranges with the same left-most and right-most positions that are found across multiple samples.
- An insert is removed if less than half of the samples have a coverage of under 20X for this insert.
- The length (in bp), the mean of each quality metric across samples, the duplication state (an insert is identified as potential duplicated or aligned multiple times in the genome based on the quality metrics) and the number of samples that contain the insert is added for each one.
- Then, the coverage per sample for the insert is added a single column at the end of the line.
- A list of all inserts (with their positions, length, mean quality metrics across samples, duplication status, number of samples and coverage for each sample) is outputted in a new file for each sequencing library.

### 05. Loci (overlapping inserts) extraction
The lists of inserts (one list per library) are used as input for the loci extraction. Loci are here defined as regions that contain all potential overlapping or unique inserts. The main steps for this script are as follows:
- The left-most position of each insert across all samples is first considered as the starting position for the loci.
	- Inserts with a left-most position in a 1,000 bp window of an existing insert (i.e. existing left-most position) are considered as the same insert. 
- Inserts are checked for overlap. If the left-most and right-most positions of two inserts are the same (+/- 20 bp), they are considered to be in the same locus.
	- Unique inserts (with no overlapping inserts) are considered to be loci. 
- A list of all loci (with their positions, number of samples that contain them and sequencing library of origin) is outputted in one file. 
- Figures are generated to investigate the complex loci (3 inserts, incomplete overlaps, etc.).

### 06. Methylation status call on loci
The list of loci (one for both libraries) is used as input for the methylotyping (methylation status calling). This methylation call is based on the comparative analysis of the sequencing libraries (see Figure 1). The main steps for this script are as follows:
- Loci unique to either library are outputted to the first output file (which corresponds to loci with only band). These loci are considered as methylated regions.
- Loci shared between libraries (that can be found in at least 1 sample in each library) are outputted to the second output file (which corresponds to loci with two bands).
- For the two bands file, methylation status are called on the loci based on the comparative analysis of the libraries.
- A list of loci with two bands (with their positions, the number of samples for each methylation status called and the methylation status for each sample) is outputted in the third output file.

