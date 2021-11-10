import os, glob, math, sys, subprocess, time, threading, random, argparse, traceback
import utilities

def extract_loci(block_dict, table_path, stab_thresh, loci_dir, all_species, win_size, win_slide, encode_multiz, stdout=True):

    from utilities\
        import num_seq_col,\
               strand_col,\
               mean_z_score_col,\
               p_score_col,\
               block_col,\
               slice_idx_col,\
               locus_idx_col,\
               species_start_col

    species_end_col = species_start_col + len(all_species)

    # Descriptions of all windows in table that are assigned to a locus
    all_win_recs = open(table_path).read().split('\n')[1:] # Read out header
    all_win_recs = (x.split('\t') for x in all_win_recs if x !='')
    all_win_recs = [(x[block_col], x[strand_col], int(x[slice_idx_col]), int(x[locus_idx_col]), x[species_start_col: species_end_col]) for x in all_win_recs if x[locus_idx_col] != 'NA']
    
    # List of triples (locus name, locus alignment file (clustal format), ungapped locus (fasta format))
    loci_alignment_list = []

    # Iterate over syntenic blocks
    block_group_list = utilities.bin_list(all_win_recs, key = lambda x: x[0])
    num_paths = len(block_group_list)
    for i, block_group in enumerate(block_group_list):
        
        block = block_group[0][0].split('.')[0]
        block_path = block_dict[block]
        
        # -------------------------------------------------------------------------------
        # EDIT: Create new variable to act as name for output directory
        block_filepath = block_group[0][0].split('.')[0]
        # -------------------------------------------------------------------------------
    
        if encode_multiz:
            # Get MAF alignments of syntenic blocks
            
            # -------------------------------------------------------------------------------
            # EDIT: Removed additional '.maf' from filepath as it is not needed.
            maf_list = [x for x in open(block_path).read().split('\n') if len(x)>0 and x[0]=='s']
            # -------------------------------------------------------------------------------
            
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
        # -------------------------------------------------------------------------------
        # EDIT: changed the location of the locus output directory, this removes the .maf from the filepath
        locus_dir = os.path.join(loci_dir, block_filepath)
        # -------------------------------------------------------------------------------
        
        if not os.path.isdir(locus_dir): os.makedirs(locus_dir)
        
        # Iterate over loci
        locus_group_list = utilities.bin_list(block_group, key = lambda x: x[3])
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
                locus_seq_list = [utilities.complement(x) for x in locus_seq_list]

            # Remove columns with only gaps
            locus_seq_list = utilities.remove_redundant_gaps(locus_seq_list)

            ### Write locus ###
            locus_idx = str(locus_group[0][3])

            # Clustal format
            clustal_path = os.path.join(locus_dir, locus_idx + '.clustal')
            clustal_string = utilities.generate_clustal(locus_header_list, locus_seq_list)
            open(clustal_path, 'w').write(clustal_string)

            # Ungapped Fasta format
            ungap_fasta_path = os.path.join(locus_dir, locus_idx + '.ungap')            
            ungap_fasta_string = '\n'.join(['>' + x + '\n' + y.replace('-', '') for x,y in zip(locus_header_list, locus_seq_list)]) + '\n'
            open(ungap_fasta_path, 'w').write(ungap_fasta_string)

            locus_name = '%s%s%s' % (block, utilities.block_locus_delim, locus_idx)
            loci_alignment_list.append([locus_name, clustal_path, ungap_fasta_path])

        
    # -------------------------------------------------------------------------------
    # EDIT: Changed variable name from 'locus_alignment_list' to 'loci_alignment_list' as 
    #       it was misspelled in REAPR v1. 
    if stdout: print >>sys.stdout, '\n'.join(['\t'.join(x) for x in loci_alignment_list])
    # -------------------------------------------------------------------------------

    return loci_alignment_list
