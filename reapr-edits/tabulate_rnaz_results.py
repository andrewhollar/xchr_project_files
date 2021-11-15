"""
Reads in the .rnaz files produced by run_RNAz_screen.py and compiles
the results into one table.  The table contains one entry for every
strand of every window in every .rnaz file.  An entry is tab-delimited
with the following columns (TODO: include a short description of each
column)

See <table_header> global variable for the column names

An .rnaz file contains the result of running RNAz with a parameter
'--both-strands' on a contig that is sliced into windows with
rnazWindow.pl.  This means that an .rnaz file starts with the RNAz
results of the forward and reverse strand of the first window and
continues to that of the last window.  The RNAz results for one strand
of a window begins with a metadata section followed by the alignment
of the sequences.  See the RNAz documentation for more details on the
output of rnazWindow.pl and RNAz).
"""

import os, glob, sys, math, argparse
import utilities

# The start and end strings that separate the RNAz output of consecutive windows
# -------------------------------------------------------------------------------
# EDIT: Changed the start_border to match the updated RNAz output. This is needed to 
#       correctly parse the RNAz output.
start_border = '############################  RNAz 2.1.1  ##############################'
# -------------------------------------------------------------------------------

end_border = '######################################################################'

def get_table_header(all_species):
    table_header = ['sequences', \
                    'columns', \
                    'read_direction', \
                    'mean_pair_identity', \
                    'shannon_entropy', \
                    'gc_content', \
                    'mean_single_seq_mfe', \
                    'consensus_mfe', \
                    'energy_contribution', \
                    'cov_contribution', \
                    'combinations', \
                    'mean_z_score', \
                    'sci', \
                    'background_model', \
                    'decision_model', \
                    'SVM_decision_value', \
                    'SVM_RNA_prob', \
                    'prediction', \
                    'mean_single_seq_mfe', \
                    'variance_single_seq_mfe', \
                    'variance_z_score', \
                    'syn_block_name', \
                    'slice_index', \
                    'locus_index']
    table_header += all_species  #[x for x in open(species).read().split('\n') if x!='']
    table_header = '#' + '\t'.join(table_header) + '\n'
    return table_header

def parse_windows(rnaz_path, log_path, block, index_path, alternate_strands, all_species):
    """
    Returns a list table_record_list where the i-th element is a list
    describing the i-th window of a syntenic block.

    Input
    --------------
    @block     : name of a syntenic block
    @rnaz_path : file path of the RNAz run over the syntenic block
    @log_path  : verbose output of the rnazWindow.pl run over the
                 syntenic block
    """

    # The start and end character sequences that separate the RNAz output of consecutive windows
    
    # -------------------------------------------------------------------------------
    # EDIT: Changed the start_border to match the updated RNAz output. This is needed to 
    #       correctly parse the RNAz output.
    start_border = '############################  RNAz 2.1.1  ##############################'
    # -------------------------------------------------------------------------------

    end_border   = '######################################################################'

    # -------------------------------------------------------------------------------
    # EDIT: Added a try/except block to assess the case where an .rnaz file is not produced.
    try:
        window_results_list = open(rnaz_path).read().split(start_border)[1:]
    except IOError:
        return []
    # -------------------------------------------------------------------------------

    win_to_slice = [int(x) for x in open(index_path).read().split('\n') if x != '']
    win_idx = 0         # Keep track of the window index
    strand = 'forward'  # Keep track of the window strand

    # Iterate over windows in the RNAz file
    table_record_list = []
    for window_results in window_results_list:

        # Parse results for window line-by-line.            
        '''
        The following is a list of hand-counted indices for the line numbers and word offsets of the window's metadata. 

        Line number | word offset | example content

        2    | 12  | Sequences: 12
        3    | 10  | Columns: 120
        4    | 20  | Reading direction: forward
        5    | 25  | Mean pairwise identity:  80.38
        6    | 18  | Shannon entropy: 0.39807
        7    | 14  | G+C content: 0.62115
        8    | 27  | Mean single sequence MFE: -40.97
        9    | 16  | Consensus MFE: -11.78
        10   | 22  | Energy contribution: -11.12
        11   | 26  | Covariance contribution:  -0.66
        12   | 20  | Combinations/Pair:   1.50
        13   | 15  | Mean z-score:  -0.22
        14   | 31  | Structure conservation index:   0.29
        15   | 19  | Background model: mononucleotide
        16   | 17  | Decision model: sequence based alignment quality
        17   | 21  | SVM decision value:  -0.51
        18   | 28  | SVM RNA-class probability: 0.287696
        19   | 13  | Prediction: OTHER
        '''

        lines = window_results.split('\n')
        sequences = lines[2][12:]
        columns = lines[3][10:]
        read_direction = lines[4][20:]
        mean_pair_identity = lines[5][25:]
        shannon_entropy = lines[6][18:]
        gc_content = lines[7][14:]
        mean_single_seq_mfe = lines[8][27:]
        consensus_mfe = lines[9][16:]
        energy_contribution = lines[10][22:]
        cov_contribution = lines[11][26:]
        combinations = lines[12][20:]
        mean_z_score = lines[13][15:]
        sci = lines[14][31:]
        background_model = lines[15][19:]
        decision_model = '_'.join(lines[16][17:].split())
        SVM_decision_value = lines[17][21:]
        SVM_RNA_prob = lines[18][28:]
        prediction = lines[19][13:]

        # Make sure that the reading direction follows the expected alternation between 'forward' and 'reverse''
        assert read_direction == strand

        # Dictionary that keeps track of which species are contained in window.
        # species_present[<species_name>] = 1 if <species_name> is in window, and 0 otherwise.
        species_present = dict.fromkeys(all_species, 0)

        # Tally which species are present
        # Parse the mfe and z_score of each single sequence
        # Calculate the variance of the mfe of single sequences and the variance of the z_score of single sequences
        variance_single_seq_mfe = 0
        variance_z_score = 0
        for i in range(lines.index(end_border) + 2, len(lines) - 2, 3):
            assert lines[i][0] == '>' # Make sure that lines[i] is the start of a sequence's data
            species = lines[i][1:].split()[0] # Extract the species name
            
            # -------------------------------------------------------------------------------
            if "6way_block_00000740" in block:
                print species
            # -------------------------------------------------------------------------------
            
            if species == 'consensus': continue # Skip the consensus sequence

            # # Hack for the Encode Multiz alignment in order to remove the chromosome name from the species identifier in the RNAz output
            # if args.encode_multiz: species = species.split('.')[0]

            # if '2R_19576798_19582822.maf.rnaz' in rnaz_path and SVM_RNA_prob == '0.667406':
            #     print >>sys.stderr, species, species_present, all_species

            species_present[species] = 1
            MFE_clause = ' '.join(lines[i+2].split(' ')[1:]) 
            MFE_clause = MFE_clause.split(',')
            mfe = float(MFE_clause[0][2:])
            z_score = float(MFE_clause[1][12:])
            variance_single_seq_mfe += (mfe - float(mean_single_seq_mfe))**2
            variance_z_score += (z_score - float(mean_z_score))**2
            
        # -------------------------------------------------------------------------------
        if "6way_block_00000740" in block:
            print species_present    
        # -------------------------------------------------------------------------------
        
        variance_single_seq_mfe = round(variance_single_seq_mfe /  float(sequences), 2)
        variance_z_score = round(variance_z_score / float(sequences), 2)
        assert sum(species_present.values()) == int(sequences)

        # Get slice index of this window
        if win_to_slice == []:  slice_idx = win_idx
        else:                   slice_idx = win_to_slice[win_idx]

        locus_idx = 'NA'  # locus_idx to be determined
        
        # Compile the window's entry
        table_record = [sequences, columns, read_direction, mean_pair_identity, shannon_entropy, \
                            gc_content, mean_single_seq_mfe, consensus_mfe, energy_contribution, \
                            cov_contribution, combinations, mean_z_score, sci, background_model, \
                            decision_model, SVM_decision_value, SVM_RNA_prob, prediction, \
                            mean_single_seq_mfe, str(variance_single_seq_mfe), str(variance_z_score), \
                            block, str(slice_idx), str(locus_idx)]
        # Add information about which species are present
        table_record.extend([str(x[1]) for x in sorted(species_present.iteritems(), key=lambda x: x[0])])
        table_record_list.append(table_record)
    
        # Update the strand and window index
        if alternate_strands:
            if strand == 'forward':
                strand = 'reverse'
            else:
                strand = 'forward'

                # Both strands of the window have now been parsed, so increment the window index
                win_idx += 1
        else:
            # Just increment the window index
            win_idx += 1

    return table_record_list

def merge_windows(table_record_list, threshold, alternate_strands, all_species):
    """
    Merge the windows

    <table_record_list>[i] := a list describing the i-th window of a
    syntenic block.
    
    Two windows are merged if the following conditions are met
    1.)  Their mean mfe z-score is <= t, where t is the stability threshold
    2.)  They are oriented in the same read direction
    3.)  The difference in their slice index is <= 3
    4.)  There is a sufficient overlap of species.  If one window has
    the set of species A and the other window has B, |A inter B| / |A union B| >= 0.5

    Windows are merged in the order of their slice index.  All of the
    windows in the forward direction are first merged, followed by
    those in the reverse direction.  Merged windows are labeled
    together as a locus and assigned a 0-based locus index.

    For every window that is assigned to a locus,
    <table_record_list>[i] is mutated in the locus_idx_column to store
    the locus index, or 'NA' if no locus is assigned.
    """

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

    # Sort the windows according to their slice index
    table_record_list.sort(key=lambda x: int(x[slice_idx_col]))

    # -------------------------------------------------------------------------------
    # Filter for windows that are below the MFE z-score threshold
    # EDIT: Change the column used for filtering. Instead of using the mean_z_score_column, use the p_score_col
    # filtered_list = [x for x in table_record_list if float(x[mean_z_score_col]) <= threshold]
    filtered_list = [x for x in table_record_list if float(x[p_score_col]) >= threshold]
    # -------------------------------------------------------------------------------

    # Merge only same strands of windows
    locus_idx = 0    # Counter for the current locus index
    boundary_unmerged_list = []
    if alternate_strands: strand_possibilities = ['forward', 'reverse']
    else: strand_possibilities = ['forward']
    for strand in strand_possibilities:
        
        strand_list = [x for x in filtered_list if x[strand_col] == strand]

        if len(strand_list) == 0:  continue

        # Merge windows
        for i, table_record in enumerate(strand_list[:-1]):
            # Assign this window a locus
            table_record[locus_idx_col] = str(locus_idx)
            
            # Compute the slice index gap to the next window
            next_record = strand_list[i+1]
            dist = int(next_record[slice_idx_col]) - int(table_record[slice_idx_col])
            
            # Do not merge current window with the next window if the
            # gap is greater than 3 or if the gap is between 2 and 3,
            # inclusive, but there is not a sufficient overlap of
            # species if one window's species is the set A, and the
            # other window's species is the set B then a sufficient
            # overlap is if |A intersect B| / |A union B| >= 0.5.  Do
            # merge if the gap is less than or equal to 2 slice
            # indices
            if dist > 3:
                locus_idx += 1
            elif dist >= 2:
                species = table_record[species_start_col : species_end_col]
                next_species = next_record[species_start_col : species_end_col]
                inter_size = sum([1 for x,y in zip(species , next_species) if int(x)==1 and int(y)==1])
                union_size = sum([1 for x,y in zip(species , next_species) if int(x)==1 or int(y)==1])
                if (float(inter_size) / union_size) < 0.5:
                    locus_idx += 1
                    boundary_unmerged_list.append(table_record)
                    boundary_unmerged_list.append(next_record)
                    boundary_unmerged_list.append(['inter: %i, union: %i, inter/union = %f\n\n' % (inter_size, union_size, float(inter_size)/union_size)])
        
        # Assign locus to last window
        strand_list[-1][locus_idx_col] = str(locus_idx)
        
        locus_idx += 1

    return boundary_unmerged_list

def write_table(table_path, rnaz_paths, block_names, log_paths, index_paths, alternate_strands, merge, threshold, species):

    table = open(table_path, 'w')
    table.write(get_table_header(species))

    # Create file containing list of windows that overlap but
    # weren't merge because of insufficient set of shared species
    if merge:  boundary_unmerged_file = open(table_path + '.unmerged', 'w')

    # Get a table entry from each RNAz file
    for i, (rnaz_path, block_name, log_path, index_path) in enumerate(zip(rnaz_paths, block_names, log_paths, index_paths)):

        # Parse the windows
        table_record_list = parse_windows(rnaz_path, log_path, block_name, index_path, alternate_strands, species)

        # Merge windows that pass the stability threshold, and record
        # windows that overlap but didn't merged because of
        # insufficient species overlap
        if merge:  boundary_unmerged_list = merge_windows(table_record_list, threshold, alternate_strands, species)
        
        # Record windows into table
        table.write('\n'.join(['\t'.join(x) for x in table_record_list]))
        if len(table_record_list) != 0: table.write('\n')
        if merge: boundary_unmerged_file.write('\n'.join(['\t'.join(x) for x in boundary_unmerged_list]))

    table.close()
    if merge: boundary_unmerged_file.close()
