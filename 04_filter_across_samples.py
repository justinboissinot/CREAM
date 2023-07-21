#!/usr/bin/env python3
"""Filter loci from both libraries

Usage:
    <program> infos_folder prefix output_folder
"""

# Modules
from collections import defaultdict
from statistics import mean
import sys
import os

# Parse user input
try:
    infos_folder = sys.argv[1]
    prefix = sys.argv[2]
    output_folder = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Read files
filelist = os.listdir(infos_folder)

inserts = defaultdict(dict)

files = [x for x in filelist if x.startswith(prefix)]

samples = []

for f in files:
    sample_name = f.replace("_filtered.infos", "")
    samples.append(sample_name)

    lines = [x.strip().split() for x in open(os.path.join(infos_folder, f)).readlines() if not x.startswith("Chrom")]

    for l in lines:

        # Remove small scaffolds
        if l[0].startswith("NW") or l[0].startswith("NC_029855.1"):
            continue

        insert = l[:3]
        insert[1:] = [int(x) for x in insert[1:]]
        insert = tuple(insert)
        data = l[3: ]

        inserts[insert][sample_name] = data

samples = sorted(samples)

with open(os.path.join(output_folder, prefix + ".tsv"), "wt") as outfile:
    header = ["Chrom", "From", "To", "Length", "Cigar", "Dist", "Score", "Duplication", "NumSamples"] + samples
    outfile.write("\t".join(header) + "\n")

    for insert in sorted(inserts):

        # Remove locus if less than half of the samples have a coverage under 20 for this locus
        min_cov = 20
        prop_cov = 0.5
        coverages = [int(x[0]) for x in list(inserts[insert].values())]
        length = insert[2] - insert[1]

        cigar = round(mean([float(x[1]) for x in list(inserts[insert].values())]), 2)
        dist = round(mean([float(x[2]) for x in list(inserts[insert].values())]), 2)
        score = round(mean([float(x[3]) for x in list(inserts[insert].values())]), 2)

        num_samples = len(coverages)
        num_samples_above_min_cov = len([x for x in coverages if x >= min_cov])

        if (num_samples_above_min_cov / num_samples) < prop_cov:
            continue

        # Decide if locus is duplicated
        duplication = "single"

        if cigar > 1.2 or dist > 4 or score > 25:
            duplication = "duplicated"

        # Build output line
        line = list(insert) + [length, cigar, dist, score, duplication]
        line.append(num_samples)

        for sample in samples:

            if not sample in inserts[insert]:
                geno = "0"

            else:
                geno = inserts[insert][sample][0]
                #geno = ":".join(loci[locus][sample])

            line.append(geno)

            #print(str(locus) + "\t" + sample + "\t" + geno)

        outfile.write("\t".join([str(x) for x in line]) + "\n")
