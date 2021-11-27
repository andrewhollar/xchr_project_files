from functools import singledispatch
import os
import sys 
import argparse
from Bio import AlignIO
# from reapr-edits.utilities import WINDOW_SIZE, WINDOW_SLIDE

WINDOW_SIZE = 120
WINDOW_SLIDE = 20
REFERENCE_SPECIES = "Homo_sapiens"


# Arguments: out BED file, Second RNAz score filter, Score difference,

def main():
    parser = argparse.ArgumentParser(description="Extract the BED locations of REAPR hits.")
    parser.add_argument("-t", "--threshold", type=float, default=0.9, help="The cutoff score for the second RNAz pass (Default: 0.9).")
    parser.add_argument("-d", "--difference", type=float, default=0.0, help="The allowed difference in RNAz scores (Default: 0).")
    parser.add_argument("-o", "--output", required=True, help="The filepath of the output bed file containing the REAPR hits.")
    parser.add_argument("-i", "--input_reapr_dir", required=True, help="The filepath to the directory containing REAPR output.")

    args = parser.parse_args()
    
    reapr_summary_path = os.path.join(args.input_reapr_dir, "summary.tab")
    assert os.path.isfile(reapr_summary_path)
    assert os.stat(reapr_summary_path).st_size != 0
    alignment_blocks_dir = os.path.join(args.input_reapr_dir, "alignments", "")
    assert os.path.isdir(alignment_blocks_dir)
    
    header = True
    header_line = ""
    
    # Column 0: Alignment block
    # Column 4: RNAz score (first pass)
    # Column 6: RNAz score (second pass)
    
    target_alignment_blocks = []
    target_block_slices = []
    
    with open(reapr_summary_path, "r") as summary_table:
        for reapr_hit_line in summary_table.readlines():
            if header:
                header_line = reapr_hit_line
                header = False
            else:
                line_tokens = reapr_hit_line.split()
                rnaz_first = float(line_tokens[4])
                rnaz_second = float(line_tokens[6])
                
                # print (rnaz_first, rnaz_second)
                
                if rnaz_second > args.threshold:
                    if rnaz_second - rnaz_first >= args.difference:
                        
                        target_alignment_blocks.append(line_tokens[0])
                        target_block_slices.append((int(line_tokens[1]), int(line_tokens[2])))
                    
    
    #We are looking to extract the start and end coordinates for each reapr hit.
    NUM_SAMPLES = 100
    block_idx = 0
    
    for block, slice_idx in zip(target_alignment_blocks, target_block_slices):
        
        #For each block
        # 1. get filtered.maf 
        # 2. use the slice_indices to extract the locus region. 
        
        
        block_filtered_path = os.path.join(alignment_blocks_dir, block.split('.')[0] + ".filtered.maf")
        locus_idx = block.split('/')[1]
        
        
        locus_start_column = slice_idx[0] * WINDOW_SLIDE
        # This is true only if the locus does not appear at the end of a block with < 120 nts (i.e. WINDOW_SIZE)
        locus_end_column = slice_idx[1] * WINDOW_SLIDE + WINDOW_SIZE
        
        
        # if slice_idx[1] == slice_idx[0]:
        
        unflanked_locus_start_pos = -1
        unflanked_locus_end_pos = -1
            
        for block in AlignIO.parse(block_filtered_path, "maf"):
            for sequence in block:
                if sequence.id.split('.')[0] == REFERENCE_SPECIES:
                    unflanked_locus_start_pos = int(sequence.annotations['start']) + len(str(sequence.seq)[:locus_start_column].replace("-", ""))
                    
                    if locus_end_column > len(str(sequence.seq)):
                        locus_end_column = len(str(sequence.seq))
                    
                    unflanked_locus_end_pos = int(sequence.annotations['start']) + len(str(sequence.seq)[:locus_end_column].replace("-", ""))
                    # unflanked_locus_start_pos = int(sequence.annotations['start'])
                    # unflanked_locus_end_pos = unflanked_locus_start_pos + int(sequence.annotations['size'])

        print(unflanked_locus_start_pos,unflanked_locus_end_pos, block)
        
        
        # unflanked_locus_start_pos = -1
        # unflanked_locus_end_pos = -1
        
        
        
        # block_windows_path = os.path.join(alignment_blocks_dir, block.split('/')[0] + ".windows")
        # locus_idx = block.split('/')[1]
        
        
        
        
        # unflanked_locus_start_pos = -1
        # unflanked_locus_end_pos = -1
        # prev_start = -1
        # multi_window_slide_offset = 0
        
        # window_idx = 0
        # for window in AlignIO.parse(block_windows_path, "maf"):
        #     if window_idx >= slice_idx[0] and window_idx <= slice_idx[1]:  
        #         if slice_idx[0] == slice_idx[1]:
        #             for sequence in window:
        #                 if sequence.id.split('.')[0] == REFERENCE_SPECIES:
        #                     unflanked_locus_start_pos = int(sequence.annotations['start'])
        #                     unflanked_locus_end_pos = unflanked_locus_start_pos + int(sequence.annotations['size'])
        #         else:
        #             for sequence in window:
        #                 if sequence.id.split('.')[0] == REFERENCE_SPECIES:
        #                     if unflanked_locus_start_pos == -1:
        #                         unflanked_locus_start_pos = int(sequence.annotations['start'])
        #                         prev_start = unflanked_locus_start_pos
        #                     else:
        #                         multi_window_slide_offset += (int(sequence.annotations['start']) - prev_start)
        #                         print("New slide offset {}".format(str(multi_window_slide_offset)))
        #                         unflanked_locus_end_pos = multi_window_slide_offset + int(sequence.annotations['size'])
        #                         prev_start = int(sequence.annotations['start'])
                
        #     window_idx += 1

                                
        # print(unflanked_locus_start_pos,unflanked_locus_end_pos, block)
                                
                            

                
                
            #     for sequence in window:
            #         if sequence.id.split('.')[0] == REFERENCE_SPECIES:
            #             print(sequence.id)
            #     print(window)
        
        

        
        # print(block)
        # print(block_windows_path)
        # print(locus_idx)
        # print(slice_idx)
        
        block_idx += 1
        if block_idx > NUM_SAMPLES:
            break
                    
    
    # print(target_alignment_blocks)
    # print(len(target_alignment_blocks))
    

if __name__=='__main__':
    main()
