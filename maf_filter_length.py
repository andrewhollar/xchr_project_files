"""
Filters out all alignment blocks that contain sequences 

python maf_filter_single_species.py path_to_maf_file.maf output.maf
"""

# pylint: disable=E0401
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys

maf_in_filepath = sys.argv[1]
maf_out_filepath = open(sys.argv[2], "w")

alignment_blocks = []
# alignment_block_lengths = {}

# Loop through each of the alignment blocks, checking the number of species contained within it (i.e. its length)
for msa in AlignIO.parse(maf_in_filepath, "maf"):

    filtered_msa_sequences = []

    # length_flag = False

    for sequence in msa:
        if int(sequence.annotations['size']) >= 50:
            filtered_msa_sequences.append(sequence)
            # length_flag = True

    # if not length_flag:
    if len(filtered_msa_sequences) >= 2:
        alignment_blocks.append(MultipleSeqAlignment(filtered_msa_sequences))
    # alignment_blocks.append(msa)

# Write all alignment blocks that have more than 1 species contained.
AlignIO.write(alignment_blocks, maf_out_filepath, "maf")
maf_out_filepath.close()

# print(maf_filepath, num_msas)
# print(sorted(alignment_block_lengths.items(), reverse=True))