#!/usr/bin/env python3
"""From 2 tables for MspI and HpaII, extract loci and their fragments

Usage:
    <program> input_msp input_hpa output_file
"""

# Modules
import matplotlib.pyplot as plt

from collections import defaultdict
import pandas as pd
import sys
import os

# Functions
def overlaps(r1, r2):
    """Return if range 1 (r1 = [from<int>, to<int>]) overlaps range 2 (r2 = [from<int>, to<int>])
    """
    return (r1[1] >= r2[0]+20) & (r2[1] >= r1[0]+20)

def get_locus_range(locus):
    """Return leftmost and rightmost nucleotide positions for a list of fragments
    """
    return [min([x[1] for x in locus]), max([x[2] for x in locus])]

# Parse user input
try:
    input_msp = sys.argv[1]
    input_hpa = sys.argv[2]
    output_file = sys.argv[3]
except:
    print(__doc__)
    sys.exit(1)

# Load data with Pandas
msp_full = pd.read_csv(input_msp, sep="\t")
hpa_full = pd.read_csv(input_hpa, sep="\t")

msp = [list(msp_full.loc[i, ["Chrom", "From", "To", "NumSamples"]]) + ["msp"] for i in range(len(msp_full))]
hpa = [list(hpa_full.loc[i, ["Chrom", "From", "To", "NumSamples"]]) + ["hpa"] for i in range(len(hpa_full))]

# Create a dictionary by chromosome and starting position
ranges = defaultdict(lambda: defaultdict(list))

for r in msp + hpa:
    chrom, fr, to, num, enz = r
    ranges[chrom][fr//1000].append(r)

locus_num = 0

with open(output_file, "wt") as outfile:
    for chrom in sorted(ranges):
        treated_bins = set()

        for _bin in sorted(ranges[chrom]):
            if _bin in treated_bins:
                continue

            data = ranges[chrom][_bin]
            next_bin = _bin + 1
            next_data = ranges[chrom][next_bin]

            # Add next consecutive bins if they contain data
            while next_data:
                treated_bins.add(next_bin)
                data += next_data
                next_bin += 1
                next_data = ranges[chrom][next_bin]

            data = sorted(data)

            while data:
                locus_num += 1
                locus_name = "locus_" + str(locus_num)
                locus = [data.pop(-1)]
                locus_range = locus[0][1: 3]

                #if not data:
                #    outfile.write("\t".join([locus_name] + [str(x) for x in locus]) + "\n")

                for i in range(len(data))[::-1]:
                    other_range = data[i][1: 3]

                    if overlaps(locus_range, other_range):
                        locus.append(data.pop(i))
                        locus_range = get_locus_range(locus)

                # Write fragments of the locus to file
                for l in sorted(locus):
                    outfile.write("\t".join([locus_name] + [str(x) for x in l]) + "\n")

                # Create figures for complex loci
                
                # TODO Genotype all the things!!!
                # Filter: Only one fragment
                if len(locus) == 1:
                    continue

                # Filter: Two fragments with identical ranges
                if len(locus) == 2:
                    if (locus[0][1] == locus[1][1]) and (locus[0][2] == locus[1][2]):
                        if (locus[0][4] != locus[1][4]):
                            #print(f"{locus_name} found one")
                            l = sorted(locus, key=lambda x: x[4])
                        continue

                # Filter: Four fragments with two pairs that have identical ranges
                if len(locus) == 4:
                    if (locus[0][1] == locus[1][1]) and (locus[0][2] == locus[1][2]):
                        if (locus[2][1] == locus[3][1]) and (locus[2][2] == locus[3][2]):
                            continue
                        
                figure_name = os.path.join("05_locus_figures", locus_name + ".png")

                fig = plt.figure(figsize=(3, 2))
                plt.plot([0, locus_range[1] - locus_range[0]], [0, 10], color="white")

                # Create figure of fragments
                y = 0

                for l in sorted(locus):
                    y += 1
                    col = "black" if l[-1] == "msp" else "red"
                    plt.plot([l[1] - locus_range[0], l[2] - locus_range[0]], [y, y], color=col, linestyle="-")
                    plt.text(0, 9.5, locus_name)
                    plt.text(0, 8.5, " ".join([str(locus[0][0]), str(locus_range[0])]))

                #plt.axis("off")
                plt.savefig(figure_name)
                plt.close(fig)

