"""
Count alignment blocks by length (i.e. the number of species in each)
"""

# pylint: disable=E0401
from Bio import AlignIO
import sys

# maf_filepath = '../output/01_preprocessing/filtered_sample_ssmaxid_chr19.maf'
maf_filepath = sys.argv[1]

num_msas = 0
bins = {0: 0,
        1: 0,
        2: 0,
        3: 0,
        4: 0,
        5: 0,
        6: 0,
        7: 0,
        8: 0,
        9: 0}
bins.setdefault(0)

for msa in AlignIO.parse(maf_filepath, "maf"):
    num_msas += 1
    
    bins[len(msa)] += 1

print(maf_filepath, bins, num_msas)