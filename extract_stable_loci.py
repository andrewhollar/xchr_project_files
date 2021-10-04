import os
import glob
import math
import sys 
import subprocess
import time
import threading
import random
import argparse
import traceback
# import utilities

def extract_loci(block_dict, table_path, stab_thresh, loci_dir, all_species, win_size, win_slide, encode_multiz, stdout=True):

    # from utilities\
    #     import num_seq_col,\
    #            strand_col,\
    #            mean_z_score_col,\
    #            p_score_col,\
    #            block_col,\
    #            slice_idx_col,\
    #            locus_idx_col,\
    #            species_start_col

    # Parsing RNAz files
    num_seq_col = 0
    strand_col = 2
    mean_z_score_col = 11
    p_score_col = 16
    block_col = 21
    slice_idx_col = 22
    locus_idx_col = 23
    species_start_col = 24
    block_locus_delim = '/'

    species_end_col = species_start_col + len(all_species)

    # Descriptions of all windows in table that are assigned to a locus
    all_win_recs = open(table_path).read().split('\n')[1:] # Read out header
    all_win_recs = (x.split('\t') for x in all_win_recs if x !='')
    all_win_recs = [(x[block_col], x[strand_col], int(x[slice_idx_col]), int(x[locus_idx_col]), x[species_start_col: species_end_col]) for x in all_win_recs if x[locus_idx_col] != 'NA']
    
    # List of triples (locus name, locus alignment file (clustal format), ungapped locus (fasta format))
    loci_alignment_list = []

    # Iterate over syntenic blocks
    block_group_list = bin_list(all_win_recs, key = lambda x: x[0])
    num_paths = len(block_group_list)
    for i, block_group in enumerate(block_group_list):
        
        block = block_group[0][0]
        try:
            block_path = block_dict[block]
        except KeyError:
            print("{} has already been removed.".format(block))
            continue

        if encode_multiz:
            # Get MAF alignments of syntenic blocks
            maf_list = [x for x in open(block_path + '.maf').read().split('\n') if len(x)>0 and x[0]=='s']
            seq_list = [x.split()[6] for x in maf_list]

            # Remove the chromosome names from the MAF headers
            # This is needed for concordance with the species guide tree for LocARNA
            header_list = [x.split()[1].split('.')[0] for x in maf_list]

        else:
            # Parse MAF alignment
            maf_list = [[a.split('\t') for a in x.split('\n') if a!='' and a[0]=='s'] for x in open(block_path).read().split('a score') if x!='']
            assert len(maf_list) == 1
            maf_list = maf_list[0]
            
            seq_list = [a[6] for a in maf_list]
            header_list = [a[1] for a in maf_list]

        # Make directory for syntenic block's loci
        locus_dir = os.path.join(loci_dir, block)
        if not os.path.isdir(locus_dir): os.makedirs(locus_dir)
        
        # Iterate over loci
        locus_group_list = bin_list(block_group, key = lambda x: x[3])
        for locus_group in locus_group_list:

            # Take the union of the species present in this locus's window
            species_present = [False for x in all_species]
            for window in locus_group:
                species_present = [bool(int(x)) or y for x,y in zip(window[4], species_present)]

            # Extract the locus's sequences and fasta headers 
            locus_header_list = []
            locus_seq_list = []
            start_slice_idx = min([window[2] for window in locus_group])
            end_slice_idx = max([window[2] for window in locus_group])
            start_column = start_slice_idx * win_slide
            end_column = end_slice_idx * win_slide + win_size
            for k, header in enumerate(header_list):
                if species_present[all_species.index(header)]:
                    locus_header_list.append(header)
                    locus_seq_list.append(seq_list[k][start_column : end_column])

            # Check to take the complement strand
            if locus_group[0][1] == 'reverse':
                locus_seq_list = [complement(x) for x in locus_seq_list]

            # Remove columns with only gaps
            locus_seq_list = remove_redundant_gaps(locus_seq_list)

            ### Write locus ###
            locus_idx = str(locus_group[0][3])

            # Clustal format
            clustal_path = os.path.join(locus_dir, locus_idx + '.clustal')
            clustal_string = generate_clustal(locus_header_list, locus_seq_list)
            open(clustal_path, 'w').write(clustal_string)

            # Ungapped Fasta format
            ungap_fasta_path = os.path.join(locus_dir, locus_idx + '.ungap')            
            ungap_fasta_string = '\n'.join(['>' + x + '\n' + y.replace('-', '') for x,y in zip(locus_header_list, locus_seq_list)]) + '\n'
            open(ungap_fasta_path, 'w').write(ungap_fasta_string)

            locus_name = '%s%s%s' % (block, block_locus_delim, locus_idx)
            loci_alignment_list.append([locus_name, clustal_path, ungap_fasta_path])

    # if stdout: print >>sys.stdout, '\n'.join(['\t'.join(x) for x in locus_alignment_list])

    return loci_alignment_list

# Taken from REAPR-utilities
def bin_list(lst, key):
    """
    Takes in a list and returns a list of bins such that every element
    x and y in a bin have key(x) == key(y).
    """
    lst.sort(key = key)
    curr_key = None
    bin_lst = []
    bin = []
    for x in lst:
        if key(x) != curr_key:
            bin = []
            curr_key = key(x)
            bin_lst.append(bin)
        bin.append(x)

    return bin_lst

#Taken from REAPR-utilities
def complement(seq):
    """
    Returns the complement strand of <seq>
    """
    
    reverse_dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', \
                    'a':'t', 't':'a', 'c':'g', 'g':'c', \
                    '-':'-', 'N':'N', 'n': 'n'}
    complement = list(seq)
    complement.reverse()
    complement = ''.join([reverse_dict[x] for x in complement])
    return complement

#Taken from REAPR-utilities
def remove_redundant_gaps(seq_list):
    """
    Remove columns of only gaps    
    """

    seq_list = [list(x) for x in seq_list]

    pop_offset = 0
    num_seqs = len(seq_list)
    col_to_del = [i for i in range(len(seq_list[0])) if sum([seq_list[j][i] !='-' for j in range(num_seqs)]) == 0 ]
    col_to_del = [i - a for i,a in zip(col_to_del, range(len(col_to_del)))]
    for i in col_to_del:
        for j in range(num_seqs):
            del seq_list[j][i]
    
    seq_list = [''.join(x) for x in seq_list]

    return seq_list

#Taken from REAPR-utilities
def generate_clustal(header_list, seq_list):
    """
    Returns an alignment in CLUSTAL format as a string.  Also removes
    any columns with only gaps.

    Input
    -----------------
    @header_list   : List of sequence headers
    @seq_list      : List of gapped sequences
    """

    num_seqs = len(seq_list)

    # Remove columns that only contain gaps
    pop_offset = 0
    col_to_del = [i for i in range(len(seq_list[0])) if sum([seq_list[j][i] !='-' for j in range(num_seqs)]) == 0 ]
    col_to_del = [i - a for i,a in zip(col_to_del, range(len(col_to_del)))]
    for i in col_to_del:
        for j in range(num_seqs):
            del seq_list[j][i]

    clustal_list = [x + '\t' + ''.join(y) for x,y in zip(header_list, seq_list)]
    clustal_string = clustal_standard_header + '\n'.join(clustal_list) + '\n'

    return clustal_string