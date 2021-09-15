"""
Filters out all sequence entries that are comprised entirely of gap characters (i.e. '-'). 
These sequences caused issues with REAPR in the LocARNA phase. If an alignment block contains such a 
sequence it is only kept if >=2 species remain after the full-gap sequence has been removed.

python maf_filter_full_gap_sequences.py path_to_maf_file.maf output.maf
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

    alignment_block = []

    for sequence in msa:
        ungapped = sequence.seq.replace("-", "").strip()

        if ungapped != "":
            alignment_block.append(sequence)
        else:
            print(sequence.id, sequence.seq, ungapped)

    if len(alignment_block) > 1:
        alignment_blocks.append(MultipleSeqAlignment(alignment_block))
    else:
        print("only 1 sequence remains in the alignment block, skipping")
    # alignment_blocks.append(msa)

# Write all alignment blocks that have more than 1 species contained.
AlignIO.write(alignment_blocks, maf_out_filepath, "maf")
maf_out_filepath.close()

# print(maf_filepath, num_msas)
# print(sorted(alignment_block_lengths.items(), reverse=True))