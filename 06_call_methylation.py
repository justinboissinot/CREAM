#!/usr/bin/env python3
"""From 2 tables for MspI and HpaII, extract loci and their fragments

Usage:
    <program> input_loci input_msp input_hpa output_file1 output_file2 output_file3
"""

# Modules
import matplotlib.pyplot as plt

from collections import defaultdict
import pandas as pd
import sys
import os

# Parse user input
try:
    input_loci = sys.argv[1]
    input_msp = sys.argv[2]
    input_hpa = sys.argv[3]
    output_file1 = sys.argv[4]
    output_file2 = sys.argv[5]
    output_file3 = sys.argv[6]
except:
    print(__doc__)
    sys.exit(1)

# Load data with Pandas
msp_full = pd.read_csv(input_msp, sep="\t")
hpa_full = pd.read_csv(input_hpa, sep="\t")


# Normalize data
msp_full_norm = msp_full.copy()
hpa_full_norm = hpa_full.copy()

for dataset in [msp_full_norm, hpa_full_norm]:
    sample_cols = list(dataset.columns)[9:]

    for col in sample_cols:
        dataset[col] = 100000 * dataset[col] / dataset[col].sum()

#print(msp_full)
#print(msp_full_norm)
#print(hpa_full)
#print(hpa_full_norm)

# Load loci
loci = [x.strip().split("\t") for x in open(input_loci).readlines()]

loci_dict = defaultdict(dict)
loci_count = defaultdict(int)

for l in loci:
    num = int(l[0].split("_")[1])
    loci_count[num] += 1
    loci_dict[l[0]][l[-1]] = l

list_samples = msp_full.columns[9:]
samples = []

for i in list_samples:
    sample_name = i.replace("CanMspI_1_", "")
    samples.append(sample_name)

with open(output_file1, "wt") as outfile1:
    header1 = ["Library", "Chrom", "From", "To", "Length", "Cigar", "Dist", "Score", "Duplication", "NumSamples"] + samples
    header2 = ["Chrom", "From", "To", "Length", "Cigar", "Dist", "Score", "Duplication", "NumSamples"] + samples
    outfile1.write("\t".join(header1) + "\n")
    with open(output_file2, "wt") as outfile2:
        outfile2.write("\t".join(header2) + "\t" + "\t".join(header2) + "\n")

        for n in sorted(loci_count):
            if loci_count[n] > 2:
                continue

            locus_name = "locus_" + str(n)

            locus = loci_dict[locus_name]

            # Cases with data in only one library
            if len(locus) == 1:
                if "msp" in locus.keys():
                    locus_infos = locus["msp"]
                    _id = locus_infos[1:4]

                    locus_cov = list(msp_full_norm.loc[(msp_full_norm["Chrom"] == _id[0])
                            & (msp_full_norm["From"] == int(_id[1]))
                            & (msp_full_norm["To"] == int(_id[2])),].iloc[0])

                    # Write outputs
                    outfile1.write("\t".join(["msp"] + [str(x) for x in locus_cov]) + "\n")

                elif "hpa" in locus.keys():
                    locus_infos = locus["hpa"]
                    _id = locus_infos[1:4]

                    locus_cov = list(hpa_full_norm.loc[(hpa_full_norm["Chrom"] == _id[0])
                            &(hpa_full_norm["From"] == int(_id[1]))
                            &(hpa_full_norm["To"] == int(_id[2])),].iloc[0])

                    # Write outputs
                    outfile1.write("\t".join(["hpa"] + [str(x) for x in locus_cov]) + "\n")

            if len(locus) == 2:

                locus_msp = locus["msp"]
                msp_id = locus_msp[1:4]

                msp_cov = list(msp_full_norm.loc[(msp_full_norm["Chrom"] == msp_id[0])
                        & (msp_full_norm["From"] == int(msp_id[1]))
                        & (msp_full_norm["To"] == int(msp_id[2])),].iloc[0])

                locus_hpa = locus["hpa"]
                hpa_id = locus_hpa[1:4]

                hpa_cov = list(hpa_full_norm.loc[(hpa_full_norm["Chrom"] == hpa_id[0])
                        & (hpa_full_norm["From"] == int(hpa_id[1]))
                        & (hpa_full_norm["To"] == int(hpa_id[2])),].iloc[0])

                # Write outputs
                outfile2.write("\t".join([str(x) for x in msp_cov + hpa_cov]) + "\n")

# Load the data from the two bands file
two_bands = pd.read_csv(output_file2, sep='\t')

# Split the dataframe into the two datasets
two_bands1 = two_bands.iloc[:, :77]
two_bands2 = two_bands.iloc[:, 77:]

# Define a function to calculate the correspondence between the coverage values
def calculate_correspondence(row):
    corr = []
    for i in range(68):
        if row[i+9] == 0 and row[i+86] == 0:
            corr.append(0)
        elif row[i+9] > 0 and row[i+86] > 0:
            corr.append(1)
        elif row[i+9] > 0 and row[i+86] == 0:
            corr.append(2)
        elif row[i+9] == 0 and row[i+86] > 0:
            corr.append(3)
        else:
            raise ValueError('Unexpected value in the coverage column')
    return corr

# Apply the function to each row of the dataframe
correspondence = two_bands.apply(calculate_correspondence, axis=1)

# Count the number of 0, 1, 2 and 3 in each row
count_0 = correspondence.apply(lambda x: x.count(0))
count_1 = correspondence.apply(lambda x: x.count(1))
count_2 = correspondence.apply(lambda x: x.count(2))
count_3 = correspondence.apply(lambda x: x.count(3))

# Concatenate the results into a new dataframe
result = pd.concat([two_bands1.iloc[:, :3], pd.DataFrame({'Num_A0-0': count_0, 'Num_A1-1': count_1, 'Num_B1-0': count_2, 'Num_D0-1': count_3}), correspondence.apply(pd.Series)], axis=1)

# Save the result to a file
result.to_csv(output_file3, sep='\t', index=False, header=list(two_bands1.columns[:3]) + ['Num_A0-0', 'Num_A1-1', 'Num_B1-0', 'Num_D0-1'] + list(two_bands1.columns[-68:]))
