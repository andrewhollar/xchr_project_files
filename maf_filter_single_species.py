"""
Filters out all alignment blocks that contain a single species (i.e. the reference species).
These alignment blocks hold no alignment information.

python maf_filter_single_species.py path_to_maf_file.maf output.maf
"""

# pylint: disable=E0401
from Bio import AlignIO
import sys

maf_in_filepath = sys.argv[1]
maf_out_filepath = open(sys.argv[2], "w")

alignment_blocks = []
# alignment_block_lengths = {}

# Loop through each of the alignment blocks, checking the number of species contained within it (i.e. its length)
for msa in AlignIO.parse(maf_in_filepath, "maf"):
    if len(msa) == 1:
        continue
    alignment_blocks.append(msa)

# Write all alignment blocks that have more than 1 species contained.
AlignIO.write(alignment_blocks, maf_out_filepath, "maf")
maf_out_filepath.close()

# print(maf_filepath, num_msas)
# print(sorted(alignment_block_lengths.items(), reverse=True))