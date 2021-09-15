"""
Count alignment blocks by length (i.e. the number of species in each)
"""

# pylint: disable=E0401
from Bio import AlignIO
import sys

# maf_filepath = '../output/01_preprocessing/filtered_sample_ssmaxid_chr19.maf'
maf_filepath = sys.argv[1]

num_msas = 0
bins = {}

for msa in AlignIO.parse(maf_filepath, "maf"):
    num_msas += 1
    
    if len(msa) in bins:
        bins[len(msa)] += 1
    else:
        bins[len(msa)] = 1

print(maf_filepath, bins, num_msas)