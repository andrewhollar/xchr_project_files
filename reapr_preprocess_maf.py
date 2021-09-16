"""
This script reads through a MAF file and distills it into individual MAF files corresponding
to the alignment blocks of the input file. These smaller files are written to an output directory
where they will be read by REAPR. An additional file mapping the alignment block names to their 
location is generated for use by the REAPR pipeline.


python reapr_preprocess_maf.py path_to_maf_file.maf
"""

# pylint: disable=E0401
from Bio import AlignIO
import sys
import os
import random

# Create an output directory to contain the files required as input to reapr
OUT_DIR = '/home/ahollar/reapr_001/alignments/'
if not os.path.exists(OUT_DIR):
    os.makedirs(OUT_DIR)

SAMPLE_DENOM = 10
MAX_SAMPLES = 2
SAMPLE_LENGTH = 500

# Method to pad an integer with zeros on the left, this returns a string of length num_positions.
def pad_int(input_int, num_positions):
    str_rep = str(input_int)
    while len(str_rep) < num_positions:
        str_rep = "0" + str_rep
    return str_rep

maf_filepath = sys.argv[1]
reapr_alignment_map = []
alignment_block_idx = 0
sample_size = 0
for msa in AlignIO.parse(maf_filepath, "maf"):
    if random.randint(1, SAMPLE_DENOM) == 1 and sample_size < MAX_SAMPLES and len(msa[0].seq) >= SAMPLE_LENGTH:
        # This should contain two pieces of information:
        #   1. the name of the alignment block
        #   2. the location of the new alignment block file (this should be in the 'alignments/' sub-directory)
        reapr_alignment_entry = []

        # create an output file name for this alignment block
        alignment_block_name = "6way_block_" + pad_int(alignment_block_idx, 8) + ".maf"
        reapr_alignment_entry.append(alignment_block_name)

        # use the file name to create a location for the output file
        alignment_block_output_location = os.path.join(OUT_DIR, alignment_block_name)
        reapr_alignment_entry.append(alignment_block_output_location)

        # write the alignment block to the output file
        # maf_out_filepath = open(alignment_block_output_location, "w")

        for sequence in msa:
            sequence.seq = str(sequence.seq).upper()


        with open(alignment_block_output_location, "w") as block_out:
            block_out.write("a score=0.00\n")
            for sequence in msa:
                block_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("s", sequence.id, sequence.annotations['start'], sequence.annotations['size'], sequence.annotations['strand'], sequence.annotations['srcSize'], str(sequence.seq).upper()))

        # AlignIO.write(msa, maf_out_filepath, "maf")
        # maf_out_filepath.close()

        reapr_alignment_map.append(reapr_alignment_entry)
        sample_size += 1

    alignment_block_idx += 1


# write the alignment block map file to the './reapr_preprocess/' directory
with open('/home/ahollar/reapr_001/alignment_blocks.txt', 'w') as blocks_out:
    for entry in reapr_alignment_map:
        blocks_out.write("{}\t{}\n".format(entry[0], entry[1]))



