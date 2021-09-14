"""
Count alignments over x in length


python maf_count_alignments.py path_to_maf_file.maf
"""

# pylint: disable=E0401
from Bio import AlignIO
import sys

maf_filepath = sys.argv[1]

alignment_block_lengths = {}

num_msas = 0
for msa in AlignIO.parse(maf_filepath, "maf"):
    if msa[0].annotations['size'] >= 15:
        num_msas += 1

    if msa[0].annotations['size'] in alignment_block_lengths:
        alignment_block_lengths[msa[0].annotations['size']] += 1
    else:
        alignment_block_lengths[msa[0].annotations['size']] = 1

print(maf_filepath, num_msas)
# print(sorted(alignment_block_lengths.items(), reverse=True))