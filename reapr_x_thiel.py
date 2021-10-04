from genericpath import isfile
import sys
import os
import multiprocessing
import argparse
import traceback
import subprocess
import random
from Bio import AlignIO
import commands
import time
import shutil

import run_first_rnaz_screen
import tabulate_rnaz_results
import extract_stable_loci

SAMPLE_DENOM = 100
MAX_SAMPLES = 5     #sys.maxsize
SAMPLE_LENGTH = 15

STABILITY_THRESHOLD = -1  # Upper threshold on mean z score
PROCESSES = 1  # Number of process to 
WINDOW_SIZE = 120
WINDOW_SLIDE = 20

random.seed(35)

def main():
    ''' Needed arguments:
            1. Input MAF file
            2. Output folder

        "Grand-Fathered" REAPR arguments:
            1. Species file
            2. Guide tree
            3. Delta
            4. Stability threshold
            5. Processes
    '''
    
    parser = argparse.ArgumentParser(description="Predict structural ncRNAs through realignment")
    parser.add_argument('-m', '--maf-file', required=True, help='The input alignment MAF file.')
    parser.add_argument('-o', '--output-folder', default=os.getcwd(), help='Directory to write the output files (Default: present working directory).')
    parser.add_argument('-s', '--species', required=True, help="Line-separated list of species in the input MAF file.")
    parser.add_argument('-g', '--guide-tree', required=True, help="Species guide tree in newick format, without branch lengths.")
    parser.add_argument('-d', '--delta', type=int, default=[20], nargs='+', help="Space separated list of realignment deviations (Default: 20).")
    parser.add_argument('-t', '--threshold', type=float, default=STABILITY_THRESHOLD, help="All regions that have a mean MFE z-score above this thresold will be filtered out (Default: -1).")
    parser.add_argument('-p', '--processes', type=int, default=PROCESSES, help="Number of cores to use for multiprocessing (Default: 1).")
    args = parser.parse_args()

    outF, errF = sys.stdout, sys.stderr
    OUT_DIR = args.output_folder
    # print(OUT_DIR)
    parent_pid = os.getpid()

    if not os.path.exists(OUT_DIR):
        os.makedirs(OUT_DIR)
    else:
        shutil.rmtree(OUT_DIR)
        os.makedirs(OUT_DIR)

    process_maf_file(args.maf_file, OUT_DIR)

    assert os.path.isfile(os.path.join(OUT_DIR, "alignment_blocks.txt")), 'Error: {0} is not a file.'.format(os.path.join(OUT_DIR, "alignment_blocks.txt"))
    alignment_blocks = [x.split('\t') for x in open(os.path.join(OUT_DIR, "alignment_blocks.txt")).read().split('\n') if x != '']
    alignment_block_dict = dict(alignment_blocks)
    alignment_block_names, alignment_block_paths = zip(*alignment_blocks)

    #Check that each of the extracted alignment blocks was successfully written to disk
    for block in alignment_block_paths:     
        assert os.path.isfile(block), 'Error: {0} is not a file'.format(block)
    
    assert os.path.isfile(args.species), 'Error: {0} is not a file'.format(args.species)
    species = sorted([x for x in open(args.species).read().split('\n') if x != ''])

    #Run initial RNAz screen

    #Setup RNAz parameters
    no_reference = False
    structural = False
    verbose = True
    both_strands = True
    alignment_format = "MAF"
    # RNAz_OUT_DIR = os.path.join(OUT_DIR, "rnaz_1")
    RNAz_args = [(alignment, no_reference, both_strands, WINDOW_SIZE, WINDOW_SLIDE, structural, commands.RNAz, commands.rnazWindow, OUT_DIR, None, alignment_format, verbose) for alignment in alignment_block_paths]
    print(errF, 'Start: First RNAz screen', get_time())
    RNAz_job_pool = multiprocessing.Pool(processes=args.processes)
    RNAz_log_list = RNAz_job_pool.map_async(run_first_rnaz_screen.run_first_rnaz_screen_MP, RNAz_args).get(99999999)
    # print(str(errF), 'End: First RNAz screen', get_time())
    print('End: First RNAz screen', get_time(), file=errF)

    #Compile table of RNAz screen results

    for block_name, block_path in alignment_block_dict.items():
        if not os.path.isfile(block_path):
            print("Popping block {} : {}".format(block_name, block_path), file=errF)
            alignment_block_dict.pop(block_name)

    alignment_block_paths = [a for a in alignment_block_paths if os.path.isfile(a)]
    RNAz_paths = [a + '.rnaz' for a in alignment_block_paths if os.path.isfile(a + '.rnaz')]
    RNAz_log_paths = [a + '.windows.log' for a in alignment_block_paths if os.path.isfile(a +'.windows.log')]
    RNAz_index_paths = [a + '.windows.indices' for a in alignment_block_paths if os.path.isfile(a + '.windows.indices')]
    initial_table = os.path.join(OUT_DIR, 'first_rnaz_screen.table')
    alternate_strands = True
    merge = True
    tabulate_rnaz_results.write_table(initial_table, RNAz_paths, alignment_block_paths, RNAz_log_paths, RNAz_index_paths, alternate_strands, merge, args.threshold, species)


    #Extract stable loci
    loci_dir = os.path.join(args.output_folder, 'loci')
    loci_alignment_list = extract_stable_loci.extract_loci(alignment_block_dict, initial_table, args.threshold, loci_dir, species, WINDOW_SIZE, WINDOW_SLIDE, True, stdout=False)

    print(loci_alignment_list)

    realignment_tables = [os.path.join(args.output_folder, 'locarna.d_%s.tab' % d) for d in args.delta]

    for delta, realign_table in zip(args.delta, realignment_tables):
        locus_names, ref_clustals, ungapped_fastas = zip(*loci_alignment_list)


    #RNAz#
    
# Method to pad an integer with zeros on the left, this returns a string of length num_positions.
def pad_int(input_int, num_positions):
    str_rep = str(input_int)
    while len(str_rep) < num_positions:
        str_rep = "0" + str_rep
    return str_rep

# Taken from REAPR-utilities.py
def get_time():
    return '\nTime:\n' + time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()) + '\nEpoch time: ' + str(time.time()) + '\n'

def process_maf_file(path_to_maf, out_dir):

    alignments_out_dir = os.path.join(out_dir, "alignments/")

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
            alignment_block_output_location = os.path.join(out_dir, "alignments/", alignment_block_name)
            reapr_alignment_entry.append(alignment_block_output_location)

            for sequence in msa:
                sequence.seq = str(sequence.seq).upper()

            with open(alignment_block_output_location, "w") as block_out:
                # print("writing alignment")
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
    with open(os.path.join(out_dir, "alignment_blocks.txt"), "w") as blocks_out: 
        for entry in reapr_alignment_map:
            blocks_out.write("{}\t{}\n".format(entry[0], entry[1]))

if __name__=='__main__':
    main()