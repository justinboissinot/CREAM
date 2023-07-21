#!/usr/bin/env python3
"""Extract infos to assess quality of ranges from bam files

Usage:
    <program> input_bam output_file
"""

# Modules
from statistics import stdev as sd
from collections import defaultdict
import pysam
import sys

# Parsing user input
try:
    input_bam = sys.argv[1]
    output_file = sys.argv[2]
except:
    print(__doc__)
    sys.exit(1)

# Read bam
samfile = pysam.AlignmentFile(input_bam, "rb")  # 'r' for read and 'b' for BAM

# Initiating dictionaries for the ranges and their info
# A defaultdict never raises a KeyError and provides a default value for the key that does not exist
ranges = defaultdict(lambda: defaultdict(list))
range_infos = defaultdict(lambda: defaultdict(list))

count = 0

# Loops through every read (segment) in the SAM file
for segment in samfile:
    s = segment.to_string()  # Converts the line of the read to a string
    if "SA:" in s or "XA" in s:
        continue  # Ignores reads with chimeric alignments (SA) or alternative hits (XA)

    s = s.split("\t")  # Splits the string of the read into a list of elements that were tab-delimited
    s = s[:9] + s[11:]  # Excludes elements of index 9 and 10 (SEQ and QUAL fields)

    if not s[1] in ["83", "99", "147", "163"]:
        continue  # Ignores reads without flags indicating they are mapped in the proper pair

    # Names the elements of the 's' list
    seqid, flag, chrom, pos, qual, cigar, _, pos2, template, dist, cigar2, cigarrev, score, _ = s

    # Keeps only loci of expected lengths
    if abs(int(template)) > 280 or abs(int(template)) < 50:
        continue

    # Takes the mapping position of the other read in the pair if the template length is negative. This is to assure
    # that both reads of the same pair have the same left-most (starting) position. A negative template length indicates
    # a read from the reverse strand and therefore that is the right-most one of the pair.
    if int(template) < 0:
        pos = pos2

    # Splits the strings relative to quality metrics and keeps only the value (last field)
    dist = dist.split(":")[-1]
    cigar2 = cigar2.split(":")[-1]
    cigarrev = cigarrev.split(":")[-1]
    score = score.split(":")[-1]

    # Converts strings to integers
    pos = int(pos) - 1  # Subtracts 1 since the SAM mapping positions are 1-based
    qual = int(qual)
    template = int(template)
    dist = int(dist)
    score = int(score)

    # Creates an index for each range. A range is defined as a pair of reads falling between the same starting (pos) and
    # ending positions (pos+template) on a same chromosome.
    index = (chrom, pos, pos + abs(template))

    # Reduces a read 's' to its quality metrics
    s = qual, cigar, template, dist, score

    # Adds the quality metrics of each read of a pair (a pair has a unique query template name or seqid) to each pair
    # that is found in a same range (index). The 'ranges' dictionary contains all the pairs of reads (and their quality
    # metrics) for each range.
    ranges[index][seqid].append(s)

    # If a range only has two pairs
    if len(ranges[index][seqid]) == 2:

        infos = list(ranges[index][seqid][0]) + list(ranges[index][seqid][1])

        for i, value in enumerate(infos):
            range_infos[index][i].append(value)

    count += 1

    if count > 1000000000:
        break

# Write to file
with open(output_file, "wt") as outfile:
    # header = "Chrom\tFrom\tTo\tCov\tQual1\tQual2\tCigar1\tCigar2\tDist1\tDist2\tScore1\tScore2\n"
    header = "Chrom\tFrom\tTo\tCov\tCigar\tDist\tScore\n"
    outfile.write(header)

    for pos in sorted(range_infos):

        infos = range_infos[pos]
        supp = []

        # Qual
        # supp.append(round(sum(infos[0]) / len(infos[0]), 2))
        # supp.append(round(sum(infos[5]) / len(infos[5]), 2))

        # Cigar
        if len(infos[1]) > 1:
            c1 = round(sd([len(x) for x in infos[1]]), 2)
            c2 = round(sd([len(x) for x in infos[6]]), 2)
            supp.append(max([c1, c2]))
        else:
            supp += [0]

        # Dist
        if len(infos[3]) > 1:
            d1 = round(sd(infos[3]), 2)
            d2 = round(sd(infos[8]), 2)
            supp.append(max([d1, d2]))
        else:
            supp += [0]

        # Score
        if len(infos[4]) > 1:
            s1 = round(sd(infos[4]), 2)
            s2 = round(sd(infos[9]), 2)
            supp.append(max([s1, s2]))
        else:
            supp += [0]

        outfile.write("\t".join([str(x) for x in list(pos) + [len(range_infos[pos][0])] + supp]) + "\n")
