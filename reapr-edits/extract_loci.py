import os, glob, math, sys, subprocess, time, threading, random, argparse, traceback
import utilities



#def write_bed_file():
    


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
    # -------------------------------------------------------------------------------
    # Example of an entry in all_win_recs.
    # ('6way_block_00000085.maf', 'forward', 3, 0, ['0', '1', '0', '1', '0', '1'])
    # -------------------------------------------------------------------------------
    
    # List of triples (locus name, locus alignment file (clustal format), ungapped locus (fasta format))
    loci_alignment_list = []

    # Iterate over syntenic blocks
    # The bin_list utility method will create seperate lists of window records according to the MAF file they originated from.
    block_group_list = utilities.bin_list(all_win_recs, key = lambda x: x[0])
    # -------------------------------------------------------------------------------
    # EDIT: removed unused variable.
    # num_paths = len(block_group_list)
    # -------------------------------------------------------------------------------

    for i, block_group in enumerate(block_group_list):
        
        # Get the name of the alignment block (including the .maf extension) 
        # Ex: 6way_block_00000085.maf
        block = block_group[0][0]
        
        # Get the filepath to the alignment block
        block_path = block_dict[block]
        
        
        # -------------------------------------------------------------------------------
        # EDIT: Create new variable to act as name for output directory
        # block_filepath = block_group[0][0].split('.')[0]
        # -------------------------------------------------------------------------------
    
        if encode_multiz:
            # Get MAF alignments of syntenic blocks
            # -------------------------------------------------------------------------------
            # EDIT: Removed additional '.maf' from filepath as it is not needed.
            maf_list = [x for x in open(block_path).read().split('\n') if len(x)>0 and x[0]=='s']
            # -------------------------------------------------------------------------------
            seq_list = [x.split()[6] for x in maf_list]
            # print seq_list

            # Remove the chromosome names from the MAF headers
            # This is needed for concordance with the species guide tree for LocARNA
            header_list = [x.split()[1].split('.')[0] for x in maf_list]
            # print header_list

            # -------------------------------------------------------------------------------
            # Edit: Added another list to get the contigs of the sequences included in the alignment block
            contig_list = [x.split()[1].split('.')[1] for x in maf_list]
            contig_list = []
            for entry in maf_list:
                contig_list.append(".".join(entry.split()[1].split('.')[1:]))
                # entry_name = ".".join(entry_name[1:])
            # print contig_list
            # -------------------------------------------------------------------------------

            maf_start_pos_list = [int(x.split()[2]) for x in maf_list]
            maf_entry_length_list = [int(x.split()[3])  for x in maf_list]
            maf_direction_list = [x.split()[4] for x in maf_list]
            maf_contig_lengths_list = [int(x.split()[5]) for x in maf_list]


        else:
            # Parse MAF alignment
            maf_list = [[a.split('\t') for a in x.split('\n') if a!='' and a[0]=='s'] for x in open(block_path).read().split('a score') if x!='']
            assert len(maf_list) == 1
            maf_list = maf_list[0]
            
            seq_list = [a[6] for a in maf_list]
            header_list = [a[1] for a in maf_list]

        alignment_block_genomic_coordinates = utilities.get_alignment_block_sequence_lengths(maf_list)
        # print alignment_block_genomic_coordinates
        
        # Make directory for syntenic block's loci
        # -------------------------------------------------------------------------------
        # EDIT: changed the location of the locus output directory, this removes the .maf from the filepath
        locus_dir = os.path.join(loci_dir, block.split('.')[0])
        
        # -------------------------------------------------------------------------------
        if not os.path.isdir(locus_dir): os.makedirs(locus_dir)

        
        # Iterate over loci
        # This separates the windows by the previously assigned locus index.
        locus_group_list = utilities.bin_list(block_group, key = lambda x: x[3])
        for locus_group in locus_group_list:
            
            locus_idx = str(locus_group[0][3])

            # -------------------------------------------------------------------------------
            # EDIT: added subdirectory to hold the BED files pertaining to this locus
            locus_bed_dir = os.path.join(locus_dir, locus_idx + "_BED_FILES")
            if not os.path.isdir(locus_bed_dir): os.makedirs(locus_bed_dir)
            # -------------------------------------------------------------------------------

            
            # Take the union of the species present in this locus's window
            species_present = [False for x in all_species]
            for window in locus_group:
                species_present = [bool(int(x)) or y for x,y in zip(window[4], species_present)]

            # Extract the locus's sequences and fasta headers 
            locus_header_list = []
            locus_seq_list = []
            
            # These are the indices of the start and end of the windows within the locus.
            start_slice_idx = min([window[2] for window in locus_group])
            end_slice_idx = max([window[2] for window in locus_group])
            
            # These are the indices of the start and end of the sequence included within the locus
            #start_column = start_slice_idx * win_slide
            #end_column = end_slice_idx * win_slide + win_size
            
            # -------------------------------------------------------------------------------
            # Number of sliding windows spanning the unfiltered block
            #alignment_length = alignment_block_genomic_coordinates[all_species[0]][1] - alignment_block_genomic_coordinates[all_species[0]][0]
            #num_slices_in_original_maf = int(math.ceil((alignment_length - (win_size - win_slide)) / float(win_slide)))
            # -------------------------------------------------------------------------------
            # BED information
            # -------------------------------------------------------------------------------
            # EDIT: add 20nt flanking regions to loci that have those regions in the initial alignment block.
            # if start_slice_idx != 0 and (num_slices_in_original_maf - end_slice_idx) > 0:  
            #     #print "adding 20nt to both ends of this locus"
            #     start_column -= win_slide
            #     end_column += win_slide
            # # -------------------------------------------------------------------------------
            # else:
            #     print "Need to retrieve nucleotides from the Genome FASTA files. %s" % (block)

            # unflanked_lengths = []
            # bed_lengths = []

            for k, species_name in enumerate(header_list):
                if species_present[all_species.index(species_name)]:
                    
                    
                    #New method for obtaining the flanking regions:
                    #1. Get the start and end positions of the aligned sequence in the MAF file
                    #       - These will be genomic coordinates
                    #2. Use BedTools/other method to flank this entire region.
                    #       - This will guarantee that whatever the locus region is within the MAF file
                    #       - will be flanked
                    #3. Determine the start column of locus within MAF alignment block
                    #4. Determine the end column of locus within MAF alignment block
                    #5. Extract flanked sequence.                    
                    # BED information
                    
                    contig_name = contig_list[k]
                    num_start = utilities.FLANK_VALUE
                    num_end = utilities.FLANK_VALUE
                    
                    flanked_maf_start_pos = maf_start_pos_list[k] - utilities.FLANK_VALUE
                    if flanked_maf_start_pos < 0:
                        num_start += flanked_maf_start_pos
                        flanked_maf_start_pos = 0
                    flanked_maf_end_pos = maf_start_pos_list[k] + maf_entry_length_list[k] + utilities.FLANK_VALUE
                    if flanked_maf_end_pos > maf_contig_lengths_list[k]:
                        num_end += (maf_contig_lengths_list[k] - flanked_maf_end_pos)
                        flanked_maf_end_pos = maf_contig_lengths_list[k]
                
                    flanked_sequence = utilities.get_flanked_sequence(species_name, contig_name, flanked_maf_start_pos, flanked_maf_end_pos, locus_bed_dir, locus_idx, seq_list[k].strip(), maf_direction_list[k], maf_contig_lengths_list[k]).lower()
                    assert seq_list[k].replace('-', '').lower() in flanked_sequence.lower()
                    
                    maf_start_column = start_slice_idx * win_slide
                    
                    if end_slice_idx - start_slice_idx > 0:
                        assert len(seq_list[k]) > 120
                        maf_end_column = end_slice_idx * win_slide + win_size
                    else:
                        if len(seq_list[k]) < 120:
                            assert start_slice_idx == 0
                            maf_end_column = len(seq_list[k])
                        else:
                            #print species_name, contig_name, len(seq_list[k]), start_slice_idx, end_slice_idx
                            maf_end_column = end_slice_idx * win_slide + win_size
                            
                            
                            print maf_end_column, len(seq_list[k])
                            
                            assert maf_end_column <= len(seq_list[k])
                    
                    #print species_name, contig_name, maf_start_column, maf_end_column 
                    
                    unflanked_region = seq_list[k][maf_start_column:maf_end_column].replace("-", "").lower()
                    # leading_region = seq_list[k][:maf_start_column].replace("-", "").lower()
                    # trailing_region = seq_list[k][maf_end_column:].replace("-", "").lower()
                    #print len(leading_region), len(trailing_region)
                    
                    #flanked_start_col = 
                    
                    #extract_from_flank_start = len(leading_region) + num
                    start_offset = flanked_sequence.find(unflanked_region)
                    #print start_offset, num_start
                    assert start_offset - num_start >= 0
                    
                    extracted_seq = flanked_sequence[(start_offset - num_start) : start_offset + len(unflanked_region) + num_end]

                    #print extracted_seq
                    
                    #maf_end_column = 
                
                    #leading_sequence = seq_list[k][:start_column].replace("-", '')
                    # print leading_sequence
                    # leading_sequence = leading_sequence.replace("-", '')
                    
                    #trailing_sequence = seq_list[k][end_column + 1:].replace("-", '')
                    # print trailing_sequence
                    # trailing_sequence = trailing_sequence.replace("-", '')

                    # bed_start = alignment_block_genomic_coordinates[species_name][0] + len(leading_sequence)
                    # bed_end = alignment_block_genomic_coordinates[species_name][1] - len(trailing_sequence)
                    
                    # unflanked_lengths.append(len(seq_list[k][start_column : end_column]))
                    # bed_lengths.append(bed_end - bed_start)
                    
                    # unflanked_seq = seq_list[k][start_column : end_column].replace("-", "")
                    
                    # print (species_name, contig_name, bed_start, bed_end)
                    
                    # If the locus is just a single window, then I have to check how long the sequence is in the maf
                    # The previous indices assumed that the window was at least 120nt (the window size)
                    # if start_slice_idx == end_slice_idx:
                    #     continue
                    # else:
                        
                    # locus_start_offset = start_slice_idx * win_slide
                    
                    
                    
                    # if end_slice_idx - start_slice_idx == 0:
                    #     locus_length = maf_entry_length_list[k] - locus_start_offset
                    # else:
                    #     locus_length = (end_slice_idx * win_slide + len(seq_list[k][locus_start_offset:]) - locus_start_offset
                    
                    # #if end_slice_idx == start_slice_idx and end_slice_idx == 0:
                    # #    locus_length = maf_entry_length_list[k]
                    # #else:
                    # #    locus_length = (end_slice_idx * win_slide + win_size) - locus_start_offset                        
                        
                    # # This isn't what I want. Instead I need to get the start/end positions of the stable locus.
                    # maf_start_pos = maf_start_pos_list[k]
                    # maf_end_pos = maf_start_pos_list[k] + maf_entry_length_list[k]
                    
                    
                    # locus_start_pos = maf_start_pos + locus_start_offset
                    # locus_end_pos = locus_start_pos + locus_length
                    
                    # flanked_start_pos = locus_start_pos - utilities.FLANK_VALUE
                    # if flanked_start_pos < 0:
                    #     flanked_start_pos = 0
                    # flanked_end_pos = locus_end_pos + utilities.FLANK_VALUE
                    # if flanked_end_pos > maf_contig_lengths_list[k]:
                    #     flanked_end_pos = maf_contig_lengths_list[k]
                    
                    # locus_start_pos = maf_start_pos + start_column
                    # locus_end_pos = maf_start_pos + end_column
                    
                    # print (species_name, contig_name, maf_start_pos, maf_end_pos, locus_start_pos, locus_end_pos, flanked_start_pos, flanked_end_pos)
                    # print (locus_start_pos, locus_end_pos)
                    
                    #locus_start_pos
                    
                    # if maf_direction_list[k] == "-":
                    #     maf_start_pos = maf_contig_lengths_list[k] - maf_start_pos_list[k]
                    #     maf_end_pos = maf_start_pos + maf_entry_length_list[k]
                    
                    #utilities.confirm_matching_sequence(species_name, contig_name, maf_start_pos, maf_end_pos, locus_bed_dir, locus_idx, seq_list[k].replace("-", "").strip(), maf_direction_list[k], maf_contig_lengths_list[k])
                    # utilities.get_flanked_sequence(species_name, contig_name, bed_start, bed_end, locus_bed_dir, locus_idx, unflanked_seq)
                    
                    
                    locus_header_list.append(species_name)
                    #locus_seq_list.append(seq_list[k][start_column : end_column])
                    locus_seq_list.append(extracted_seq)
                    
            # print unflanked_lengths
            # print bed_lengths

            # Check to take the complement strand
            if locus_group[0][1] == 'reverse':
                locus_seq_list = [utilities.complement(x) for x in locus_seq_list]

            # Remove columns with only gaps
            # -------------------------------------------------------------------------------
            # There should not be any gaps at this point
            # locus_seq_list = utilities.remove_redundant_gaps(locus_seq_list)
            # -------------------------------------------------------------------------------

            # if "6way_block_00000740" in block:
            #     print species_present
            #     print maf_list
                #print locus_seq_list


            ### Write locus ###

            # Clustal format
            # clustal_path = os.path.join(locus_dir, locus_idx + '.clustal')
            # clustal_string = utilities.generate_clustal(locus_header_list, locus_seq_list)
            # open(clustal_path, 'w').write(clustal_string)

            # Ungapped Fasta format
            ungap_fasta_path = os.path.join(locus_dir, locus_idx + '.ungap')            
            ungap_fasta_string = '\n'.join(['>' + x + '\n' + y.replace('-', '') for x,y in zip(locus_header_list, locus_seq_list)]) + '\n'
            open(ungap_fasta_path, 'w').write(ungap_fasta_string)

            locus_name = '%s%s%s' % (block, utilities.block_locus_delim, locus_idx)
            # old : loci_alignment_list.append([locus_name, clustal_path, ungap_fasta_path])
            loci_alignment_list.append([locus_name, ungap_fasta_path])
        
    # -------------------------------------------------------------------------------
    # EDIT: Changed variable name from 'locus_alignment_list' to 'loci_alignment_list' as 
    #       it was misspelled in REAPR v1. 
    if stdout: print '\n'.join(['\t'.join(x) for x in loci_alignment_list])
    # -------------------------------------------------------------------------------

    return loci_alignment_list
