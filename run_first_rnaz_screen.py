import os
import sys
import subprocess
import math
import traceback
import time
import threading
import shutil

#Taken from REAPR-run_RNAz_screen.py
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

#Modified from REAPR-run_RNAz_screen.py
def run_rnazWindow(align_path, windows_path, verbose_path, no_reference, rnazWindow_command, window_size, window_slide, verbose=False):
    """ 
    Runs rnazWindow.pl on the alignment in <align_path> and outputs
    sliced windows to <windows_path
    """

    no_reference = '--no-reference' if no_reference else ''
    cmd = '%s --min-length=50 --window=%s --slide=%s --max-gap=0.7 --verbose --max-seqs=40 --min-id=0 %s %s' % (rnazWindow_command, str(window_size), str(window_slide), no_reference, align_path)

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

#Modified from REAPR-run_RNAz_screen.py
def run_RNAz(windows_path, rnaz_path, both_strands, structural_model, rnaz_command, verbose=False):
    """
    Runs RNAz on the sliced alignments in <windows_path> and outputs
    to <rnaz_path>.
    """

    log = ''

    both_strands = '--both-strands' if both_strands else ''
    structural_model = '-l' if structural_model else ''
    cmd = '%s -d %s %s --cutoff=0 %s' % (rnaz_command, structural_model, both_strands, windows_path)

    # Open output files with 1mb buffer
    rnaz_output = open(rnaz_path, 'w', 1000000)    
    
    # Run RNAz
    start_time = time.time()
    p = subprocess.Popen(cmd, shell=True, stdout=rnaz_output, stderr=subprocess.PIPE)
    stderr = p.communicate()[1]
    if stderr!='': log += 'RNAz stderr:\n' + str(stderr)  # RNAz's stderr

    # Print running time
    log += '\n' + cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds'
    
    rnaz_output.close()

    return log

def run_first_rnaz_screen(alignment, no_reference, both_strands, window_size, window_slide, structural, RNAz, rnazWindow, out_dir, tmp_dir, alignment_format, verbose=False):
    if tmp_dir is None:
        tmp_dir = out_dir
    
    log = ''
    name = os.path.basename(alignment)
    length = get_seq_length(alignment, form=alignment_format)

    # Number of sliding windows spanning the block
    num_slices = int(math.ceil((length - (window_size - window_slide)) / float(window_slide)))

    # Run rnazWindow
    windows_path = os.path.join(tmp_dir, "alignments/", name + '.windows')
    verbose_path = os.path.join(tmp_dir, "alignments/", name + '.windows.log')
    log += run_rnazWindow(alignment, windows_path, verbose_path, no_reference, rnazWindow, window_size, window_slide, verbose)

    # Write the window to slice index map
    win_to_slice_path = os.path.join(tmp_dir, "alignments/", name + '.windows.indices')
    index_windows(verbose_path, [], num_slices, win_to_slice_path)

    # Run RNAz
    rnaz_path = os.path.join(tmp_dir, "alignments/", name + '.rnaz')
    log += '\n' + run_RNAz(windows_path, rnaz_path, both_strands, structural, RNAz, verbose)

    if verbose: print(str(sys.stderr), log, "\n")
    return log

def run_first_rnaz_screen_MP(jobs):
    try:
        return run_first_rnaz_screen(*jobs)
    except KeyboardInterrupt:
        pass
    except:
        print(sys.stderr, traceback.print_exc())

#Taken from REAPR-utilities.py
def get_seq_length(path, form):
    curr_align = open(path).read()
    assert curr_align != '', 'Alignment is empty'

    form = form.upper()

    if form=='MAF':

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