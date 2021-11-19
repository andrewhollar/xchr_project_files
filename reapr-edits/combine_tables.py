import sys, os, argparse, bisect, subprocess
import utilities 
import commands

def run_alistat(alistat, alignment):
    """ <alignment> needs to be CLUSTAL format. """

    alistat_cmd = '%s -q %s' % (alistat, alignment)
    stdout, stderr = subprocess.Popen(alistat_cmd, shell=True, stdout=subprocess.PIPE).communicate()
    identity = [float(x.split()[2].split('%')[0]) / 100. for x in stdout.split('\n') if x != '' and x.split()[0] == 'Average' and x.split()[1] == 'identity:']
    try:
        assert len(identity) == 1
        identity = identity[0]
    except AssertionError:
        identity = 'NA'
    return identity

def run_compalignp(compalignp, ref_path, test_path):
    cmd = '%s -r %s -t %s' % (compalignp, ref_path, test_path)
    stdout = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE).communicate()[0]
    try:
        return float(stdout)
    except:
        print >>sys.stderr, stdout
        
        # -------------------------------------------------------------------------------
        # EDIT: Commented out the following line, as it appears this was left in erroneously
        #0 / asdf
        # -------------------------------------------------------------------------------

def combine_tables(original_table, realign_tables, deltas, loci_dir, guide_tree, species, alistat, compalignp, output_table):

    from utilities\
        import num_seq_col,\
               strand_col,\
               mean_z_score_col,\
               p_score_col,\
               block_col,\
               slice_idx_col,\
               locus_idx_col,\
               species_start_col
               
    # original table = 'original_wga.tab'
    # realign tables =  'locarna.d_%s.tab' where %s = args.delta
    # loci_dir = /reapr_x/loci
    # guide_tree = True
    # species = list of species
    # output_table = 'summary.tab'
    

    # Read initial table. Skip windows not assigned to a locus for not
    # passing stability threshold
    initial_lines = [x.split() for x in open(original_table).read().split('\n')[1:]
                        if x!='' and x[0]!='#' and x.split()[locus_idx_col]!='NA']
    locus_names = sorted(set(['%s%s%s' % (line[block_col], utilities.block_locus_delim, line[locus_idx_col]) for line in initial_lines]))
    print locus_names
    # raise IOError("End")


    # Make dictionary mapping locus name to a list of features:
    # -- min slice index
    # -- max slice index
    # -- read direction, 0 = forward, 1 = backward
    # -- RNAz p-score before realignment (max score in windows)
    # -- alistat score before realignment
    # for i in deltas:
        # -- RNAz p-score after realignment w/ locarna (max score in windows)
        # -- alistat score after realignment w/ locarna
        # -- compalignp score for change w/ locarna
    # -- number of sequences
    # -- sequence bit vector

    loci_dict = dict([(x, [-1, -1, 'NA', 0, 'NA'] +\
                          reduce(lambda a,b:a+b, [[0, 'NA', 'NA'] for i in range(len(realign_tables))], []) +\
                          ['NA']+\
                          [0 for a in range(len(species))])\
                      for x in locus_names])

    print loci_dict

    # Dictionary : delta --> list of lines from realignment table.
    # Also, update RNAz score after realignment
    realign_dict = dict()
    for i, (d, table) in enumerate(zip(deltas, realign_tables)):
        realign_lines = [x.split() for x in open(table).read().split('\n')[1:] if x!='' and x[0]!='#']
        for line in realign_lines:
            features = loci_dict[line[block_col]]
            features[5 + 3*i] = max(features[5 + 3*i], float(line[p_score_col]))

        realign_dict[d] = realign_lines

    # For later naming
    realign_suffix = 'locarna%s' % ('.g' if guide_tree else '')
        
    for line in initial_lines:

        locus_name = '%s%s%s' % (line[block_col], utilities.block_locus_delim, line[locus_idx_col])
        
        # Many lines aren't considered because the window wasn't
        # assigned to a locus, or the locus couldn't be realigned
        if loci_dict.has_key(locus_name): features = loci_dict[locus_name]
        else:                             continue

        # Update slice index range
        slice_idx = int(line[slice_idx_col])
        features[0] = slice_idx if (features[0]==-1) else min(features[0], slice_idx)
        features[1] = slice_idx if (features[1]==-1) else max(features[1], slice_idx)

        # Update locus read direction
        if features[2]=='NA':  features[2] = line[strand_col]
        else:                  features[2] == line[strand_col]

        # Update RNAz score before realignment
        p_score = float(line[p_score_col])
        features[3] = max(features[3], p_score)

        if alistat:
            # Update alistat score before realignment
            
            # -------------------------------------------------------------------------------
            # EDIT: changed 'locus_name' to not include the '.maf' extension
            original = os.path.join(loci_dir, locus_name.replace('.maf', '') + '.clustal')
            print original
            # -------------------------------------------------------------------------------    
        
            if features[4] == 'NA':
                features[4] = run_alistat(commands.alistat, original)

        for i, d in enumerate(deltas):
            
            # -------------------------------------------------------------------------------
            # EDIT: changed 'locus_name' to not include the '.maf' extension
            realignment = os.path.join(loci_dir, locus_name.replace('.maf', '') + '.clustal.%s.%s' % (realign_suffix,d))
            # -------------------------------------------------------------------------------


            if alistat:
                # Update alistat score after realignment
                if features[6 + 3*i]=='NA':
                    features[6 + 3*i] = run_alistat(commands.alistat, realignment)

            if compalignp:
                # Update compalignp score
                if features[7 + 3*i]=='NA':
                    features[7 + 3*i] = run_compalignp(commands.compalignp, original, realignment)

        # Update species bit vector and number of sequences
        start_col, end_col, species_end_col = 6+3*len(deltas), 6+3*len(deltas)+len(species), species_start_col+len(species)
        features[start_col: end_col] = [int(x) or y for x,y in zip(line[species_start_col: species_end_col], features[start_col: end_col])]
        features[5+3*len(deltas)] = sum(features[start_col: end_col])

    loci = [list(x) for x in loci_dict.items()]    

    # Sort the loci.  Recall that to sort by multiple features, sort
    # in order of increasingly important features (think radix sort)
    
    # Sort by locus index
    loci.sort(key=lambda a: int(a[0].split(utilities.block_locus_delim)[1]))
    # Sort by wga block
    loci.sort(key=lambda a: a[0].split(utilities.block_locus_delim)[0])
    # Sort by max difference in RNAz scores before vs after realignment
    loci.sort(key=lambda a: - max([a[1][5+3*i] - a[1][3] for i in range(len(deltas))]))

    # Convert to strings
    for a in loci: a[1] = map(str, a[1])
    
    # Write table
    table_header = ['locus_name', 'start_slice_idx', 'end_slice_idx', 'strand', 'RNAz_score.original', 'alistat.original'] +\
                   reduce(lambda a,b:a+b, [['%s.%s.%s' % (x, realign_suffix, d) for x in ['RNAz_score', 'alistat', 'compalignp']] for d in deltas], []) +\
                   ['no_of_species'] + species
    assert len(table_header) == 1 + (6 + 3*len(deltas) + len(species))
    table_header = '#' + '\t'.join(table_header)
    open(output_table, 'w').write(table_header + '\n' + '\n'.join([x[0]+'\t'+'\t'.join(x[1]) for x in loci]) + '\n')
