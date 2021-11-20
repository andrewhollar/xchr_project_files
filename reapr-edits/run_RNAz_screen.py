import subprocess, argparse, os, glob, sys, math, threading, traceback, time, random, multiprocessing, threading, shutil
import utilities

def index_windows(log_path, other_removals, num_slices, win_to_slice_path):
    """
    1.) Reads in the log output of rnazWindow.pl when run with the
    --verbose option 

    2.) Infers what are the 0-based slice indices of
    the windows outputted by rnazWindow.pl based on the list of
    discarded windows in the log output.

    3.) <other_removals> is a list of window indices removed because
    of the removal of consensus gap sequence.  See function
    remove_consensu_gap_seq.

    4.) Returns a list <win_to_slice>[i] := the slice index of the
    i-th window returned by rnazWindow.pl on a syntenic block
    """

    verbose_log = [x for x in open(log_path).read().split('\n') if x !='']
    # Output nothing if verbose log is empty
    if len(verbose_log) == 0: 
        open(win_to_slice_path, 'w').write('')
    else:
        # Compute <win_to_slice>[i] = slice index of the i-th window output by rnazWindows.pl

        # -------------------------------------------------------------------------------
        # EDIT: Attempting to update the discard_idx_list, I suspect that this no longer works due
        #       to changes made to the output of rnazWindow.pl
        #discard_idx_list = [verbose_log[i-1].split(':')[0].split()[-1] for i,x in enumerate(verbose_log) if 'discarded' in x]
        # discard_idx_list = [verbose_log[i-2].split(':')[1].split()[-1] for i,x in enumerate(verbose_log) if 'discarded' in x]
        
        discard_idx_list = []
        
        for i, line in enumerate(verbose_log):
            if 'discarded' in line:
                if 'window' in verbose_log[i - 1]:
                    discard_idx_list.append(verbose_log[i - 1].split(':')[0].split()[-1])
                elif 'discarded' in verbose_log[i - 1]:
                    continue
                else:
                    discard_idx_list.append(verbose_log[i - 2].split(':')[1].split()[-1])
                        
        # example #1
#         Outside training range: Alignment 1, window 15: Seq 3: base composition out of range.
#         Removing Seq 3
#         Alignment 1 discarded: No sequences left.
        
        # example #2
        # Alignment 1, window 45: Removing seq 2: too many gaps.
        # Alignment 1 discarded: Too few sequences left.
        
        
        
        # -------------------------------------------------------------------------------        
        
        try:
            assert len(discard_idx_list) == len(set(discard_idx_list))
            # Offset the indices by 1 to make them 0-based
            discard_idx_list = [int(x) - 1 for x in discard_idx_list]
            # # Bad implementation if verbose log doesn't have a line for the last window
            # # Total number of sliced windows in contig
            # for i in range(len(verbose_log)-1, -1, -1):
            #     if verbose_log[i].split()[2] == 'window':
            #         num_slices = int(verbose_log[i].split(':')[0].split()[3])
            #         break
            win_to_slice = sorted(list(set(range(num_slices)) - set(discard_idx_list)))

            # Remove other window indices that were removed after rnazWindows.pl.
            # For example, because of the removal of consensus gap seq
            for pop_offset, i in enumerate(other_removals):
                win_to_slice.pop(i - pop_offset)

            open(win_to_slice_path, 'w').write('\n'.join([str(x) for x in win_to_slice]) + '\n')
        except AssertionError:
            print log_path, discard_idx_list

# -------------------------------------------------------------------------------
# EDIT: Add this method which acts as a filter prior to running rnazWindow. This 
#       uses the rnazSelectSeqs command to remove sequences that have 100% MPI with
#       the reference species.
def run_rnazSelectSeqs(alignment_path, rnazSelectSeqs_output_path, rnazSelectSeqs_command):

    cmd = '%s --max-id=99 %s' % (rnazSelectSeqs_command, alignment_path)
    
    rnazSelectSeqs_output = open(rnazSelectSeqs_output_path, 'w', int(1e6)) 
    
    start_time = time.time()
    subprocess.Popen(cmd, shell=True, stdout=rnazSelectSeqs_output, stderr=subprocess.PIPE).wait()
    log = cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds\n'

    rnazSelectSeqs_output.close()

    return log
# -------------------------------------------------------------------------------

def run_rnazWindow(align_path, windows_path, verbose_path, no_reference, rnazWindow_command, window_size, window_slide, verbose=False):
    """ 
    Runs rnazWindow.pl on the alignment in <align_path> and outputs
    sliced windows to <windows_path
    """

    log = []

    no_reference = '--no-reference' if no_reference else ''
    # -------------------------------------------------------------------------------
    # EDIT: 1) Added the --max-gap=0.7 option to remove sequences with too many gap characters.
    #       2) Changed the --min-id value from 0 to 40 to remove alignment blocks with too low a 
    #       mean pairwise identity.
    #       3) Added the --min-length=50 option to remove sequences below that length.
    #       4) Added the --max-masked=1.0 option to match Clayton's script.
    #       Both changes were done to match the Thiel publication.

    cmd = '%s --min-length=50 --window=%s --slide=%s --verbose --max-seqs=40 --min-id=40 --max-gap=0.7 --max-masked=1.0 %s %s' % (rnazWindow_command, str(window_size), str(window_slide), no_reference, align_path)
    # -------------------------------------------------------------------------------

    # Open output files with 1mb buffer
    verbose_log = open(verbose_path, 'w', int(1e6)) 
    window_output = open(windows_path, 'w', int(1e6))

    # Run rnazWindow.pl
    start_time = time.time()
    subprocess.Popen(cmd, shell=True, stdout=window_output, stderr=verbose_log).wait()

    # Print running time
    #log = '\n' + cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds'
    log.append(cmd)
    log.append('Running time: ' + str(time.time() - start_time) + ' seconds')

    verbose_log.close()
    window_output.close()

    return log

def run_RNAz(windows_path, rnaz_path, both_strands, structural_model, rnaz_command, verbose=False):
    """
    Runs RNAz on the sliced alignments in <windows_path> and outputs
    to <rnaz_path>.
    """

    log = [] #''

    both_strands = '--both-strands' if both_strands else ''
    structural_model = '-l' if structural_model else ''
    #cmd = '%s -d %s %s --show-gaps --cutoff=0 %s' % (rnaz_command, structural_model, both_strands, windows_path)

    # -------------------------------------------------------------------------------
    # EDIT: 1) Removed the '--show-gaps' command line parameter as it is no longer supported in RNAz 2.1.1
    #       2) Changed the cutoff from 0.0 to 0.5 to only report those that pass the threshold. 
    #          >> This is likely to cause problems with the rest of the pipeline.
    cmd = '%s -d %s %s --cutoff=0.0 %s' % (rnaz_command, structural_model, both_strands, windows_path)
    # -------------------------------------------------------------------------------

    # Open output files with 1mb buffer
    rnaz_output = open(rnaz_path, 'w', 1000000)    
    
    # Run RNAz
    start_time = time.time()
    p = subprocess.Popen(cmd, shell=True, stdout=rnaz_output, stderr=subprocess.PIPE)
    stderr = p.communicate()[1]
    if stderr!='': log += 'RNAz stderr:\n' + stderr  # RNAz's stderr

    # Print running time
    #log += '\n' + cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds'
    log.append(cmd)
    log.append('Running time: ' + str(time.time() - start_time) + ' seconds')

    
    rnaz_output.close()

    return log

# -------------------------------------------------------------------------------
# EDIT: Added the rnazSelectSeqs argument, which provides the location of the perl script.
# -------------------------------------------------------------------------------
def eval_alignment(alignment, no_reference, both_strands, window_size, window_slide, structural, RNAz, rnazWindow, rnazSelectSeqs, out_dir, tmp_dir, alignment_format, pass_idx, verbose=False):
    """
    Run an RNAz screen on a MAF alignment
       1.) Slice alignment into windows with rnazWindows.pl
       2.) Run RNAz on window

    Input
    -------------------------------------------
    alignment       : File path to MAF alignment
    log             : File handler to logging file
    window_size     : The number of columns in each window
    window_slide    : The number of columns to move
    RNAz            :
    rnazWindow      :
    out_dir         :
    verbose         :
    no_reference    : 

    Output
    ------------------------------------------
    out_dir/<alignment>.windows
    out_dir/<alignment>.windows.log
    out_dir/<alignment>.windows.indices
    out_dir/<alignment>.rnaz
    """

    if out_dir is None:  out_dir = os.path.dirname(alignment)

    if tmp_dir is None:
        redirect = False
        tmp_dir = out_dir        
    else:
        redirect = True
        
        # Make another temporary subdirectory to avoid collisions over
        # same locus numbers (but different blocks) during multiprocessing
        x = utilities.get_random_string(8)
        tmp_dir = os.path.join(tmp_dir, x)
        os.mkdir(tmp_dir)

    # Initialize the log variable to hold the output from rnazWindow and RNAz procedures.
    log = [] #''

    alignment_name = os.path.basename(alignment)

    # Length of alignment (i.e. # of columns)
    alignment_length = utilities.get_seq_length(alignment, form=alignment_format)

    # Number of sliding windows spanning the block
    num_slices = int(math.ceil((alignment_length - (window_size - window_slide)) / float(window_slide)))

    if pass_idx == 1:
        # -------------------------------------------------------------------------------
        # EDIT: Run rnazSelectSeqs
        filtered_maf_path = os.path.join(tmp_dir, alignment_name[:-4] + '.filtered.maf')
        log.append(run_rnazSelectSeqs(alignment, filtered_maf_path, rnazSelectSeqs))
        # -------------------------------------------------------------------------------   

        # Run rnazWindow
        windows_path = os.path.join(tmp_dir, alignment_name + '.windows')
        verbose_path = os.path.join(tmp_dir, alignment_name + '.windows.log')
        # -------------------------------------------------------------------------------
        # EDIT: Changed the input alignment to be the output from rnazSelectSeqs.
        #if pass_idx == 1:
        log.extend(run_rnazWindow(filtered_maf_path, windows_path, verbose_path, no_reference, rnazWindow, window_size, window_slide, verbose))
        # else:
        #     log += run_rnazWindow(alignment, windows_path, verbose_path, no_reference, rnazWindow, window_size, window_slide, verbose)
        # -------------------------------------------------------------------------------

        # Write the window to slice index map
        win_to_slice_path = os.path.join(tmp_dir, alignment_name + '.windows.indices')
        index_windows(verbose_path, [], num_slices, win_to_slice_path)

        # Run RNAz
        rnaz_path = os.path.join(tmp_dir, alignment_name + '.rnaz')
        
        # -------------------------------------------------------------------------------
        # EDIT: Only run RNAz if there has been information extracted about the windows.
        if not os.stat(windows_path).st_size == 0:
            log.extend(run_RNAz(windows_path, rnaz_path, both_strands, structural, RNAz, verbose))
            #log += '\n' + run_RNAz(windows_path, rnaz_path, both_strands, structural, RNAz, verbose)
        # -------------------------------------------------------------------------------

        if redirect:
            for suffix in ['.windows', '.windows.log', '.rnaz', '.windows.indices']:
                dest = os.path.join(out_dir, alignment_name + suffix)
                if os.path.isfile(dest): os.remove(dest)   # remove before writing
                shutil.move(os.path.join(tmp_dir, alignment_name + suffix), dest)
            shutil.rmtree(tmp_dir)

        # if verbose: print >>sys.stderr, log + '\n'
    elif pass_idx == 2:
        # Do not need to run RNAzSelectSeqs, or RNAzwindow.
        
        # Run RNAz
        rnaz_path = os.path.join(tmp_dir, alignment_name + '.rnaz')
        # -------------------------------------------------------------------------------
        # EDIT: Only run RNAz if there has been information extracted about the windows.
        # if not os.stat(windows_path).st_size == 0:
        
        # print alignment_length
        
        if alignment_length < 400 and alignment_length > 49:
            log.append("running rnaz on realigned locus of length %s" % (str(alignment_length)))
            log.extend(run_RNAz(alignment, rnaz_path, both_strands, structural, RNAz, verbose))
                    
        if alignment_length >= 400:
            log.append("rnaz skipped because the improved boundaries of locus excede a width of 400 nt (%s)" % (alignment_length))
        elif alignment_length <= 49:
            log.append("rnaz skipped because the improved boundaries of locus failed to excede a width of 50 nt (%s)" % (alignment_length))

    #if verbose: print log + '\n'

    return log


# This method is the main driver of the first RNAz screen. This will begin the procedure on
# each of the alignment blocks of the MAF MSA.
def eval_alignment_multiprocessing(x):
    try:
        return eval_alignment(*x)
    except KeyboardInterrupt:
        pass
    except:
        print traceback.print_exc()
