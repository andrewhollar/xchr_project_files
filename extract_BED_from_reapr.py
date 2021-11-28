import os
import sys 
import argparse
from Bio import AlignIO

WINDOW_SIZE = 120
WINDOW_SLIDE = 20
FLANK_VALUE = 20
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
    loci_dir = os.path.join(args.input_reapr_dir, "loci", "")
    assert os.path.isdir(alignment_blocks_dir)
    
    header = True
    
    # Column 0: Alignment block
    # Column 4: RNAz score (first pass)
    # Column 6: RNAz score (second pass)
    
    target_alignment_blocks = []
    target_block_slices = []
    
    with open(reapr_summary_path, "r") as summary_table:
        for reapr_hit_line in summary_table.readlines():
            if header:
                header = False
            else:
                line_tokens = reapr_hit_line.split()
                rnaz_first = float(line_tokens[4])
                rnaz_second = float(line_tokens[6])                
                if rnaz_second > args.threshold:
                    if rnaz_second - rnaz_first >= args.difference:
                        target_alignment_blocks.append(line_tokens[0])
                        target_block_slices.append((int(line_tokens[1]), int(line_tokens[2])))
                    
    #We are looking to extract the start and end coordinates for each reapr hit.
    NUM_SAMPLES = 100
    block_idx = 0
    
    BED_COORDINATES = []
    
    for block, slice_idx in zip(target_alignment_blocks, target_block_slices):
        block_filtered_path = os.path.join(alignment_blocks_dir, block.split('.')[0] + ".filtered.maf")
        locus_idx = block.split('/')[1]
        locus_improved_boundaries_path = os.path.join(loci_dir, block.split('.')[0], locus_idx + ".ungap.locarna.g.20.d", "reliability_fit_once.out")        
        
        reliability_profile_fit_once_output = open(locus_improved_boundaries_path).read().split('\n')
        fit_line = ""
        for line in reliability_profile_fit_once_output:
            if "FIT" in line:
                fit_line = line    
        fit_line_tokens = fit_line.split()
        assert fit_line_tokens[0] == "FIT"
        fit_line_tokens = fit_line_tokens[1:]
        
        fit_line_indices = [int(x) for x in fit_line_tokens]
        boundaries_start = min(fit_line_indices)
        boundaries_end = max(fit_line_indices)
                
        locus_start_column = slice_idx[0] * WINDOW_SLIDE
        # This is true only if the locus does not appear at the end of a block with < 120 nts (i.e. WINDOW_SIZE)
        locus_end_column = slice_idx[1] * WINDOW_SLIDE + WINDOW_SIZE

        
        unflanked_locus_start_pos = -1
        unflanked_locus_end_pos = -1
        
        flanked_locus_start_pos = -1
        flanked_locus_end_pos = -1

        contig_name = ""
        seq_strand = ""
            
        for maf_block in AlignIO.parse(block_filtered_path, "maf"):
            for sequence in maf_block:
                if sequence.id.split('.')[0] == REFERENCE_SPECIES:
                    contig_name = sequence.id.split('.')[1]              
                    if sequence.annotations['strand'] == -1:
                        seq_strand = "-"
                    else:
                        seq_strand = "+"
                    
                    unflanked_locus_start_pos = int(sequence.annotations['start']) + len(str(sequence.seq)[:locus_start_column].replace("-", ""))
                    if locus_end_column > len(str(sequence.seq)):
                        locus_end_column = len(str(sequence.seq))
                    unflanked_locus_end_pos = int(sequence.annotations['start']) + len(str(sequence.seq)[:locus_end_column].replace("-", ""))
                    
                    flanked_locus_start_pos = unflanked_locus_start_pos - FLANK_VALUE
                    if flanked_locus_start_pos < 0:
                        flanked_locus_start_pos = 0
                    flanked_locus_end_pos = unflanked_locus_end_pos + FLANK_VALUE
                    if flanked_locus_end_pos > int(sequence.annotations['srcSize']):
                        flanked_locus_end_pos = int(sequence.annotations['srcSize'])

        improved_boundaries_start = flanked_locus_start_pos + boundaries_start
        improved_boundaries_end = flanked_locus_start_pos + boundaries_end

        bed_entry = "\t".join([contig_name, str(improved_boundaries_start), str(improved_boundaries_end), block.split('.')[0], "0", seq_strand])
        BED_COORDINATES.append(bed_entry)
        
        block_idx += 1
        # if block_idx > NUM_SAMPLES:
        #     break
        
    
    with open(args.output, "w") as bed_out:
        for entry in BED_COORDINATES:
            bed_out.write(entry + "\n")
    

if __name__=='__main__':
    main()
