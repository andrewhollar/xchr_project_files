import os
import random
import sys
import shutil
from Bio import AlignIO

SAMPLE_DENOM = 100
MAX_SAMPLES = 2    #sys.maxsize
SAMPLE_LENGTH = 15
random.seed(35)

# -------------------------------------------------------------------------------
# EDIT: Added the following function to generate the required files from the input MAF file
def main(path_to_maf, OUT_DIR):
    
    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        shutil.rmtree(OUT_DIR)
        os.makedirs(OUT_DIR)

    
    alignments_out_dir = os.path.join(OUT_DIR, "alignments/")

    if not os.path.exists(alignments_out_dir):
        os.makedirs(alignments_out_dir)

    reapr_alignment_map = []
    alignment_block_idx = 0
    sample_size = 0
    for msa in AlignIO.parse(path_to_maf, "maf"):
        if random.randint(1, SAMPLE_DENOM) == 1 and sample_size <= MAX_SAMPLES and len(msa[0].seq) >= SAMPLE_LENGTH:
            # This should contain two pieces of information:
            #   1. the name of the alignment block
            #   2. the location of the new alignment block file (this should be in the 'alignments/' sub-directory)
            reapr_alignment_entry = []

            # create an output file name for this alignment block
            alignment_block_name = "6way_block_" + pad_int(alignment_block_idx, 8) + ".maf"
            reapr_alignment_entry.append(alignment_block_name)

            # use the file name to create a location for the output file
            alignment_block_output_location = os.path.join(OUT_DIR, "alignments/", alignment_block_name)
            reapr_alignment_entry.append(alignment_block_output_location)

            for sequence in msa:
                sequence.seq = str(sequence.seq).upper()

            with open(alignment_block_output_location, "w") as block_out:
                block_out.write("a\tscore=0.00\n")
                for sequence in msa:

                    if sequence.annotations['strand'] == 1:
                        strand_char = "+"
                    else:
                        strand_char = "-"

                    block_out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format("s", sequence.id, sequence.annotations['start'], sequence.annotations['size'], strand_char, sequence.annotations['srcSize'], str(sequence.seq).upper()))

            reapr_alignment_map.append(reapr_alignment_entry)
            sample_size += 1

        if sample_size > MAX_SAMPLES:
            break

        alignment_block_idx += 1

    # write the alignment block map file to the './reapr_preprocess/' directory
    with open(os.path.join(OUT_DIR, "alignment_blocks.txt"), "w") as blocks_out: 
        for entry in reapr_alignment_map:
            blocks_out.write("{}\t{}\n".format(entry[0], entry[1]))
# -------------------------------------------------------------------------------

# Method to pad an integer with zeros on the left, this returns a string of length num_positions.
def pad_int(input_int, num_positions):
    str_rep = str(input_int)
    while len(str_rep) < num_positions:
        str_rep = "0" + str_rep
    return str_rep

if __name__=='__main__':
    main(sys.argv[1], sys.argv[2])