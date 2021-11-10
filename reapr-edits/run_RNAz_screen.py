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
        discard_idx_list = [verbose_log[i-1].split(':')[0].split()[-1] for i,x in enumerate(verbose_log) if 'discarded' in x]
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


def run_rnazWindow(align_path, windows_path, verbose_path, no_reference, rnazWindow_command, window_size, window_slide, verbose=False):
    """ 
    Runs rnazWindow.pl on the alignment in <align_path> and outputs
    sliced windows to <windows_path
    """

    no_reference = '--no-reference' if no_reference else ''
    cmd = '%s --window=%s --slide=%s --verbose --max-seqs=40 --min-id=0 %s %s' % (rnazWindow_command, str(window_size), str(window_slide), no_reference, align_path)

    # Open output files with 1mb buffer
    verbose_log = open(verbose_path, 'w', int(1e6)) 
    window_output = open(windows_path, 'w', int(1e6))

    # Run rnazWindow.pl
    start_time = time.time()
    subprocess.Popen(cmd, shell=True, stdout=window_output, stderr=verbose_log).wait()

    # Print running time
    log = cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds'

    verbose_log.close()
    window_output.close()

    return log

def run_RNAz(windows_path, rnaz_path, both_strands, structural_model, rnaz_command, verbose=False):
    """
    Runs RNAz on the sliced alignments in <windows_path> and outputs
    to <rnaz_path>.
    """

    log = ''

    both_strands = '--both-strands' if both_strands else ''
    structural_model = '-l' if structural_model else ''
    cmd = '%s -d %s %s --show-gaps --cutoff=0 %s' % (rnaz_command, structural_model, both_strands, windows_path)

    # Open output files with 1mb buffer
    rnaz_output = open(rnaz_path, 'w', 1000000)    
    
    # Run RNAz
    start_time = time.time()
    p = subprocess.Popen(cmd, shell=True, stdout=rnaz_output, stderr=subprocess.PIPE)
    stderr = p.communicate()[1]
    if stderr!='': log += 'RNAz stderr:\n' + stderr  # RNAz's stderr

    # Print running time
    log += '\n' + cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds'
    
    rnaz_output.close()

    return log

def eval_alignment(alignment, no_reference, both_strands, window_size, window_slide, structural, RNAz, rnazWindow, out_dir, tmp_dir, alignment_format, verbose=False):
    """
    Run an RNAz screen on a MAF alignment
       1.) Slice alignment into windows with rnazWindows.pl
       2.) Run RNAz on window

    Input
    -------------------------------------------
    alignment       : File path to MAF alignment
    log             : File handler to logging file
    window_size     :
    window_slide    :
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

    log = ''

    name = os.path.basename(alignment)

    # Length of alignment (i.e. # of columns)
    length = utilities.get_seq_length(alignment, form=alignment_format)

    # Number of sliding windows spanning the block
    num_slices = int(math.ceil((length - (window_size - window_slide)) / float(window_slide)))

    # Run rnazWindow
    windows_path = os.path.join(tmp_dir, name + '.windows')
    verbose_path = os.path.join(tmp_dir, name + '.windows.log')
    log += run_rnazWindow(alignment, windows_path, verbose_path, no_reference, rnazWindow, window_size, window_slide, verbose)

    # Write the window to slice index map
    win_to_slice_path = os.path.join(tmp_dir, name + '.windows.indices')
    index_windows(verbose_path, [], num_slices, win_to_slice_path)

    # Run RNAz
    rnaz_path = os.path.join(tmp_dir, name + '.rnaz')
    log += '\n' + run_RNAz(windows_path, rnaz_path, both_strands, structural, RNAz, verbose)

    if redirect:
        for suffix in ['.windows', '.windows.log', '.rnaz', '.windows.indices']:
            dest = os.path.join(out_dir, name + suffix)
            if os.path.isfile(dest): os.remove(dest)   # remove before writing
            shutil.move(os.path.join(tmp_dir, name + suffix), dest)
        shutil.rmtree(tmp_dir)

    if verbose: print >>sys.stderr, log + '\n'
    return log

def eval_alignment_multiprocessing(x):
    try:
        return eval_alignment(*x)
    except KeyboardInterrupt:
        pass
    except:
        print >>sys.stderr, traceback.print_exc()
