import time, subprocess, os, glob, math, pprint, random, shutil, bisect, sys, traceback, string

def get_time():
    return '\nTime:\n' + time.strftime("%a, %d %b %Y %H:%M:%S", time.gmtime()) + '\nEpoch time: ' + str(time.time()) + '\n'

def fasta_to_clustal(fasta_paths, clustal_paths):
    """
    Converts FASTA alignments into CLUSTAL format.  Also removes columns with only gaps.

    <fasta_paths> := list of the paths of FASTA alignments
    <clustal_paths> := list of the corresponding output paths to write CLUSTAL alignments
    """

    for fasta_path, clustal_path in zip(fasta_paths, clustal_paths):
        fasta_list = [x.split('\n') for x in open(fasta_path).read().split('>') if x != '']
        header_list = [x[0] for x in fasta_list]
        seq_list = [list(''.join(x[1:])) for x in fasta_list]
        num_seqs = len(seq_list)

        # Remove columns that only contain gaps
        pop_offset = 0
        col_to_del = [i for i in range(len(seq_list[0])) if sum([seq_list[j][i] !='-' for j in range(num_seqs)]) == 0 ]
        col_to_del = [i - a for i,a in zip(col_to_del, range(len(col_to_del)))]
        for i in col_to_del:
            for j in range(num_seqs):
                del seq_list[j][i]

        clustal_list = [x + '\t' + ''.join(y) for x,y in zip(header_list, seq_list)]
        clustal_string = '\n'.join(clustal_list)

        # Write clustal header
        clustal_file = open(clustal_path, 'w')       
        clustal_file.write(clustal_standard_header)
        clustal_file.write(clustal_string)
        clustal_file.write('\n')
        clustal_file.close()

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

def clustal_to_fasta(clustal_paths, fasta_paths):
    """
    Converts CLUSTAL alignments into FASTA format

    <clustal_paths> := list of the paths of CLUSTAL alignments
    <fasta_paths> := list of the corresponding output paths to write FASTA alignments
    """

    # Change into lists
    if not isinstance(clustal_paths, type([])):
        clustal_paths = [clustal_paths]
    if not isinstance(fasta_paths, type([])):
        fasta_paths = [fasta_paths]
    
    for clustal_path, fasta_path in zip(clustal_paths, fasta_paths):
        header_order = []
        header_to_seq = dict()
        for line in open(clustal_path).read().split('\n'):
            if line == '': continue
            elif line[0] == '#': continue
            elif line[:7] == 'CLUSTAL': continue
            else:
                header, seq = line.split()
                if header not in header_order: header_order.append(header)
                if header_to_seq.has_key(header):
                    header_to_seq[header].append(seq)
                else:
                    header_to_seq[header] = [seq]
        f = open(fasta_path, 'w')
        for header in header_order:
            seq = ''.join(header_to_seq[header])
            f.write('>%s\n%s\n' % (header, seq))
        f.close()

def clustal_to_fasta_string(clustal):
    """
    Same as clustal_to_fasta, but does conversion for a single clustal
    alignment as a string object, rather than file.  Returns a string
    of the FASTA alignment.
    """
    header_order = []
    header_to_seq = dict()
    for line in clustal.split('\n'):
        if line == '': continue
        elif line[0] == '#': continue
        elif line[:7] == 'CLUSTAL': continue
        else:
            header, seq = line.split()
            if header not in header_order: header_order.append(header)
            if header_to_seq.has_key(header):
                header_to_seq[header].append(seq)
            else:
                header_to_seq[header] = [seq]
    fasta_list = []
    for header in header_order:
        seq = ''.join(header_to_seq[header])
        fasta_list.append('>%s\n%s\n' % (header, seq))
    return ''.join(fasta_list)

def fix_clustal_header(clustal_path):
    '''
    Replaces the header line of a clustal alignment from a nonstandard
    variation like LocARNA's clustal header to the standard header
    found in constants.py.
    '''

    clustal_lines = open(clustal_path).read().split('\n')
    assert clustal_lines[0][:len('CLUSTAL')] == 'CLUSTAL'
    clustal_lines[0] = clustal_standard_header
    open(clustal_path, 'w').write('\n'.join(clustal_lines))

def fasta_to_maf(fasta_path, output_path):
    """
    Converts a FASTA alignment into MAF format as outlined by UCSC.
    Writes only the MAF features that are read by RNAz and ignores the
    rest.

    @param fasta_path File containing FASTA alignment.
    @param output_path File to output MAF alignment.
    """

    fasta_list = open(fasta_path).read().split('>')[1:]
    output = open(output_path, 'w')
    output.write('a score=0\n')
    for fasta in fasta_list:
        # UCSC's specs for MAF format
        # http://genome.ucsc.edu/FAQ/FAQformat.html#format5
        src = fasta.split('\n')[0]
        text = ''.join(fasta.split('\n')[1:])
        size = str(len([x for x in text if x != '-']))
        start = '0'
        srcSize = size
        strand = '+'
        line = '\t'.join(['s',src,start,size,strand,srcSize,text,'\n'])
        output.write(line)
    output.close()
      
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

def get_seq_length(path, form):
    """Calculates the length of every sequence (including gaps) in an
    alignment.    Checks that the length of all sequence in the file are
    the same and returns this uniform length.
    
    @ path   : file path to alignment
    @ form : format of alignment (either FASTA, CLUSTAL, or MAF)"""

    curr_align = open(path).read()
    assert curr_align != '', 'Alignment is empty'

    form = form.upper()
    if form=='FASTA':

        seq_lengths = [sum([len(y) for y in x.split('\n')[1:]]) for \
                           x in curr_align.split('>') if x!='']

        # Ensure equal lengths
        assert seq_lengths.count(seq_lengths[0]) == len(seq_lengths)

        seq_length = int(seq_lengths[0])

    elif form=='CLUSTAL':

        clustal_list = [x.split() for x in curr_align.split('\n') \
                            if x!='' and x[:7]!='CLUSTAL']
        clustal_list = [x for x in clustal_list if len(x) != 0]
        header_list, seq_list = [], []
        for x in clustal_list:
            if x[0] in header_list:
                seq_list[header_list.index(x[0])].append(x[1])
            else:
                header_list.append(x[0])
                seq_list.append([x[1]])
        seq_list = [''.join(x) for x in seq_list]
        seq_lengths = [len(x) for x in seq_list]

        # Ensure equal lengths
        assert seq_lengths.count(seq_lengths[0]) == len(seq_lengths)

        seq_length = int(seq_lengths[0])

    elif form=='MAF':

        # List of all alignments in file.  An alignment is a list of gapped sequences.
        maf_list = [[a.split('\t')[6] for a in x.split('\n') if a!='' and a[0]=='s'] for x in curr_align.split('a score') if x!='']
        
        # Length of every sequence in every alignment
        maf_lengths = [[len(a) for a in x] for x in maf_list]
        
        # Check that all sequences in the same alignment have same length
        assert all([all([a==x[0] for a in x]) for x in maf_lengths])
                   
        # Sum up alignment lengths
        seq_length = sum([x[0] for x in maf_lengths])
        
    else:
        raise Exception('Format not yet supported')

    return seq_length

def is_int(x):
    """
    Returns where object <x> represents an integer.  For example, the
    string "locarna" would return false, but the string "10" would
    return true.
    """
    try:
        int(x)
        return True
    except:
        return False

def partition_list(lst, parts):
    """
    Partitions the list <lst> into <parts> sub_lists of almost equal
    length.  The first <parts> - 1 sublists will be of equal length,
    and the last sublist will take the reminder.  For example, if
    lst=[1,2,3,4,5] and parts = 3, then [[1,2],[3,4],[5]] is returned.
    """

    # The length of the first parts -1 sublists
    section = int(math.ceil(len(lst) / float(parts)))
    
    partitions = []
    for i in range(parts):
        partitions.append(lst[section*i:min(len(lst), section*(i+1))])
    return partitions

def get_entropy(window, dna, normalize = True):
    """
    Returns the entropy of the window.
   
    normalize is a boolean for whether the entropy is divided by the
    number of columns.
    """

    # Do an import to the local namespace for efficiency
    from math import log
    
    # Only deal with DNA, otherwise would not need to count U
    assert dna
    # Make sure this is not RNA
    try:
        assert 'U' not in window[0]
        assert 'u' not in window[0]
    except:
        print 'window'
        print window
        assert False

    # Uppercase everything for counting
    seq_list = [seq.upper() for seq in window]
    
    num_seq = len(window)
    num_col = len(window[0])
    letters = ['A', 'T', 'C', 'G', '-', 'N']
    
    # Change window from a list of sequences to a list of alignment column
    col_list = [[seq[i] for seq in window] for i in range(num_col)]
    count_list = [[float(x.count(l)) for l in letters] for x in col_list]
    # Distribute the N count as 1/4 to every nucleotide
    count_list = [[a + .25 * x[5] for a in x[:5]] for x in count_list]
    # Get probabilities
    prob_list = [[a / num_seq for a in x] for x in count_list]
    
    # Calculate entropy
    entropy = sum([sum([0 if x == 0 else (-1 * x * log(x,2)) for x in y]) for y in prob_list]) 
    if normalize: entropy = entropy / float(num_col)

    return entropy

def run_SISSIz(window, log_file, signature, flanks):
    """
    Runs SISSIz on a CLUSTAL alignment inputted as a string.
    """
    
    start_time = time.time()
    
    # Write window to a file
    clustal_path = os.path.join(ram_dir, 'tmp%i' % random.randint(0,10000000000))
    open(clustal_path, 'w').write(window)

    threshold = 720

    while (flanks <= threshold):
        cmd = '%s -s -d --dna --clustal --flanks=%i %s' % (sissiz_command, flanks, clustal_path)
        x = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr = subprocess.PIPE)
        shuffle_window, stderr = x.communicate()
        if 'ERROR: Negative logarithm while generating distance matrix.' in stderr:
            flanks = flanks * 2
        else:  break

    # Print out message if flanks had to be raised
    if flanks > 180:
#        print 'Error with', signature, ' ; flanks', flanks
#        print cmd
        pass

    # No shuffling was possible, so just use original window
    if 'ERROR: Negative logarithm while generating distance matrix.' in stderr:
        log_file.write('No shuffling possible at %s\n' % signature)
        log_file.write('Reason: not enough flanking\n')
#        print 'No shuffling possible', signature
#        print 'Reason: not enough flanking'
        shuffle_window = window

    # Remove tmp file
#    os.remove(clustal_path)

    log_file.write('%s\n%s\nRunning time: %f seconds\n\n' % (cmd, signature, time.time() - start_time))

    return shuffle_window

def run_multiperm(window, log_file, signature):
    """
    Runs multiperm on a CLUSTAL alignment inputted as a string.
    """
    
    # Write window to a file
    clustal_path = os.path.join(ram_dir, 'tmp%i' % random.randint(0,10000000000))
    open(clustal_path, 'w').write(window)

    # Need to change directory because that's where multiperm writes
    cmd = 'cd %s ; %s -w -c level1 %s' % (ram_dir, multiperm_command, clustal_path)
    start_time = time.time()
    subprocess.Popen(cmd, shell=True).wait()
    output_path = os.path.join(ram_dir, 'perm_001_%s' % os.path.basename(clustal_path))
    shuffle_window = open(output_path).read()

    # Delete files
    os.remove(output_path)
    os.remove(clustal_path)

    log_file.write('%s\n%s\nRunning time: %f seconds\n\n' % (cmd, signature, time.time() - start_time))

    return shuffle_window    

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

def open_fasta_file(fasta_path):
    return open_fasta_string(open(fasta_path).read())

def open_fasta_string(fasta_string):
    fasta_list = [x.split('\n') for x in fasta_string.split('>') if x!='']
    header_list = [x[0] for x in fasta_list]
    seq_list = [''.join(x[1:]) for x in fasta_list]
    return header_list, seq_list

def condense_ali_comp_table(table_dir):
    paths = [x for x in glob.glob(os.path.join(table_dir, '*')) if os.path.basename(x).isdigit()]
    condense_file = open(os.path.join(table_dir, 'table'), 'w')
    for path in paths:
        condense_file.write(open(path).read())
        condense_file.write('\n')
    condense_file.close()

def read_table_with_header(f,header='all',linefilter=None):
    """ read table with column header
    and return as list of hashs.
    Project to header (which is useful to save space)
    """
    
    # print "Read table "+filename
    
    # if f is a string, try to open the file f
    filename_given=isinstance(f,str)
    if filename_given:
        f = open(f,"r")

    lines=[]

    x=f.readline()
    file_header=x.split()
    if header=='all':
        header=file_header
        
    # NOTE: don't use functional syntax here to enable space saving
    # when reducing to a smaller number of columns
    for x in f:
        if (x != '') and (linefilter==None or linefilter(x.split())):
            line_hash = dict(zip(file_header,x.split()))
            line_hash_proj =  dict([(k,line_hash[k]) for k in header])
            lines.append(line_hash_proj)

    if filename_given: f.close()

    return (header,lines)

def write_table_with_header(f,table,header='all'):
    """ write table with given header to file of given name or open file handle
    """
    
    filename_given=isinstance(f,str)
    if filename_given:
        f = open(f,"w")

    if header=='all': header=table.keys()

    # write header
    f.write('\t'.join([col for col in header])+'\n')
    
    # write table rows
    for line in table:
        f.write('\t'.join([line[col] for col in header])+'\n')

    if filename_given: f.close()
          
def col_to_pos(block, block_start, start_col, end_col, pos_to_col_dict, return_NA=False):
    """Converts a range of alignment column coordinates, start_col and
    end_col, in the Fly PECAN Alignments to a range of genomic
    coordinates in D. melanogaster.  A start_col that has a gap in
    D. mel is converted to the first genomic position after the gap.
    An end_col that has a gap is converted to the first position
    before the gap. start_col and end_col are 0-based half-closed."""

    pos_to_col = pos_to_col_dict[block]

    i = bisect.bisect_left(pos_to_col, start_col)
    j = bisect.bisect_right(pos_to_col, end_col)

    # No match to D.mel
    if i==len(pos_to_col) or j==0:
        if return_NA:
            return ('NA', 'NA')
        else:
            assert False

    try:
        start_pos = block_start + (i if pos_to_col[i]==start_col else i+1)
        end_pos = block_start + (j if (len(pos_to_col)!=j and pos_to_col[j]==end_col) else j-1)    
    except:
        traceback.print_exc()
        print 'params', block, block_start, start_col, end_col
        print 'dmel pos to col goes from', pos_to_col[0], pos_to_col[-1]
        print i, 'length of pos to col', len(pos_to_col)

    return (start_pos, end_pos)

def map_pos_to_cols(block_names, wga_dir, species, encode_multiz, verbose=False):
    """
    Maps coordinates of the genome with FASTA header <species> to the
    alignment columns in the original wga.  For the Fly WGA, you will
    probably want <species> to be 'DroMel_CAF1', and for the encode WGA, 'human'.
    
    @input
    block_names := list of the names of syntenic blocks in original wga
    
    @returns
    pos_to_col_dict := dictionary where
    
    pos_to_col_dict[block] := pos_to_col
    block := name of a syntenic block
    pos_to_col[i] := 0-based alignment column corresponding to 0-based position i of D. mel in the block
    """

    if verbose: print 'Starting mapping alignment columns to positions...\nTotal blocks: %i\nOn block number: ' % len(block_names)

    # Convert D. mel coordinates to alignment columns.
    pos_to_col_dict = dict()
    for k, block in enumerate(block_names):

        # Get sequence of reference species
        if encode_multiz:
            ref_seq = [list(x.split()[6]) for x in open(os.path.join(wga_dir, block)).read().split('\n') if len(x)>0 and x[:2+len(species)]=='s %s' % species]
        else:
            ref_seq = [list(''.join(x.split('\n')[1:])) for x \
                        in open(os.path.join(wga_dir, block)).read().split('>') \
                        if x!='' and x.split('\n')[0] == species]

        if len(ref_seq)==0:
            continue
        else:
            assert len(ref_seq) == 1
        ref_seq = ref_seq[0]

        # Update dictionary
        pos_to_col_dict[block] = [i for i, c in enumerate(ref_seq) if c != '-']
        
        # Print progress after iterating over every 1% of blocks 
        if verbose and (k % (len(block_names) / min(100, len(block_names)))) == 0: 
            print k,
            sys.stdout.flush()

    if verbose: print '   Finished mapping'

    return pos_to_col_dict

def get_species(species_path, wga_dir, form):
    """Get the name of species in a WGA from a precomputed file.  If
    the file doesn't exist, then create it by reading the WGA.
    Returns a list of the species represented in a whole genome
    alignment (WGA)
    
    @wga_dir   : Directory of whole genome alignment
    @form    : Format of alignment.  Currently supports only FASTA
    """
    if os.path.isfile(species_path):
        species = dict([x.split() for x in open(species_path).read().split('\n') if x!=''])
    else:
        form = form.upper()
        assert form in ['FASTA'], 'Support currently only for FASTA alignments'
        if form=='FASTA':
            p = subprocess.Popen("cat %s | grep '>'" % os.path.join(wga_dir, '*.fa'),shell=True, stdout=subprocess.PIPE)
            species = [x[1:] for x in set(p.communicate()[0])]
            assert p.returncode==0
        open(species_path, 'w').write('\n'.join(species))

    return species

def get_alignment_length(alignment, form):
    pass

def get_random_string(N):
    return ''.join(random.choice(string.ascii_uppercase + string.digits) for x in range(N))

def kill_group_processes():
    p = subprocess.Popen('kill $(pgrep -g %s)' % (os.getpid()), shell=True, stdout=open(os.devnull,'w'), stderr=subprocess.STDOUT)
    p.wait()
    assert p.returncode == 0

def file_exists(f):
    assert os.path.isfile(f), 'Error: {0} is not a file'.format(f)

def directory_exists(f):
    assert os.path.isdir(f), 'Error: {0} is not a directory'.format(f)


# -------------------------------------------------------------------------------
# EDIT: Added method to retrieve the genomic coordinates of the sequences included in the alignment block.
def get_alignment_block_sequence_lengths(maf_sequence_lines):
    species_genomic_coordinates = {}
    
    for maf_sequence in maf_sequence_lines:
        maf_sequence_line_tokens = maf_sequence.split('\t')
        species_name = maf_sequence_line_tokens[1].split('.')[0]
        species_start = int(maf_sequence_line_tokens[2])
        species_end = species_start + int(maf_sequence_line_tokens[3])
        
        species_genomic_coordinates[species_name] = (species_start, species_end)
    
    return species_genomic_coordinates

# EDIT: Added method to retrieve the flanked sequence of nucleotides.
def get_flanked_sequence(species, contig, start, end, locus_bed_dir, locus_idx):
    from commands import BEDTOOLS
    # print species_to_genome_dict['Homo_sapiens']
    
    bed_filepath = os.path.join(locus_bed_dir, locus_idx + "." + species + ".bed")
    bed_entry = "\t".join([contig, str(start), str(end)])
    open(bed_filepath, "w").write(bed_entry)
    
    flanked_output = os.path.join(locus_bed_dir, locus_idx + "." + species + ".flanked.fa")
    
    cmd = '%s getfasta -fi %s -fo %s -bed %s' % (BEDTOOLS, species_to_genome_dict[species], flanked_output, bed_filepath)
    start_time = time.time()
    subprocess.Popen(cmd, shell=True, stdout=subprocess.subprocess.STDOUT, stderr=subprocess.PIPE).wait()
    print 'Running time: ' + str(time.time() - start_time) + ' seconds'

# -------------------------------------------------------------------------------


GENOME_DIR = '/home/ahollar/alignments/7way_genomes'
species_to_genome_dict = {"Homo_sapiens": os.path.join(GENOME_DIR, "GCA_000001405.27_GRCh38.p12_genomic.fna"), \
                        "Macaca_mulatta": os.path.join(GENOME_DIR, "GCF_000772875.2_Mmul_8.0.1_genomic.fna"), \
                        "Callithrix_jacchus": os.path.join(GENOME_DIR, "GCA_002754865.1_ASM275486v1_genomic.fna"), \
                        "Mus_musculus": os.path.join(GENOME_DIR, "GCF_000001635.26_GRCm38.p6_genomic.fna"), \
                        "Canis_lupus_familiaris": os.path.join(GENOME_DIR, "GCF_000002285.3_CanFam3.1_genomic.fna"), \
                        "Sus_scrofa": os.path.join(GENOME_DIR, "GCF_000003025.5_Sscrofa10.2_genomic.fna"), \
                        "Oryctolagus_cuniculus": os.path.join(GENOME_DIR, "GCF_000003625.3_OryCun2.0_genomic.fna")}


# Fixed parameters: do not change #
#----------------------------------

# Parsing RNAz files
num_seq_col = 0
strand_col = 2
mean_z_score_col = 11
p_score_col = 16
block_col = 21
slice_idx_col = 22
locus_idx_col = 23
species_start_col = 24

# A standard header for CLUSTAL alignments
clustal_standard_header = 'CLUSTAL 2.0.10 multiple sequence alignment\n\n'

# Parameters for sliding window across WGA
WINDOW_SIZE = 120   # Length of Window

# -------------------------------------------------------------------------------
# EDIT: Changed the window-slide value to 20 to match the Thiel publication.
WINDOW_SLIDE = 20   # Length of slide from window to window
# -------------------------------------------------------------------------------
        
STABILITY_THRESHOLD = -1  # Upper threshold on mean z score

PROCESSES = 1  # Number of process to parallelize

# Delimiter character to join wga block names and loci numbers to form loci names
block_locus_delim = '/'

#RNAZ_WINDOW = '/home/mikeyu/bio/RNAz/share/RNAz/perl/rnazWindow2.pl'
