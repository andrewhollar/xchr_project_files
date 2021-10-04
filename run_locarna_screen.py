import os
import sys
import subprocess
import time
import traceback
import shutil
# import utilities

def run_locarna(locarna, clustal_path, ungap_fasta_path, target_dir, target_file, max_diff, alifold_consensus_dp, guide_tree, verbose=True):
    """Run LocARNA.  Returns True if successful, or False otherwise"""

    if not os.path.isdir(target_dir): os.makedirs(target_dir)

    # Setup parameters
    acd_arg          = '--alifold-consensus-dp' if alifold_consensus_dp else ''
    guide_tree_arg   = '--treefile {0}'.format(guide_tree) if guide_tree else ''
    max_diff_arg     = '--max-diff {0}'.format(max_diff)
    ref_arg          = '--max-diff-aln={0}'.format(clustal_path)
    target_dir_arg   = '--tgtdir={0}'.format(target_dir)

    # Run LocARNA, (don't include --write-structure because it's unnecessary since intermediate directory will be deleted?)
    cmd = '%s --probabilistic --consistency_transformation %s %s %s %s --keep-sequence-order %s %s' \
        % (locarna, acd_arg, guide_tree_arg, ref_arg, max_diff_arg, target_dir_arg, ungap_fasta_path)

    start_time = time.time()
    subprocess.Popen(cmd, shell=True, stdout=open('/dev/null','w'), stderr=subprocess.STDOUT).wait()
    if verbose: print(cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds\n', file=sys.stderr)

    # Copy final alignment 'results/result.aln' and delete
    # intermediate directory unless final alignment wasn't created,
    # probably due to impossible banding constraints
    result_path = os.path.join(target_dir, 'results', 'result.aln')
    if os.path.isfile(result_path):
        
        fix_clustal_header(result_path)  # Fix clustal header (remove non-standard LocARNA header)
        shutil.copyfile(result_path, target_file)   # Copy final alignment
        shutil.rmtree(target_dir)                   # Delete intermediate directory

        return True
    else:
        if verbose: print('Error: no alignment could be produced.', file=sys.stderr)
        return False

def run_locarna_pool(x):
    try:
        return run_locarna(*x)
    except KeyboardInterrupt:
        pass
    except:
        print(traceback.print_exc(), file=sys.stderr)

#Taken from REAPR-utilities
def fix_clustal_header(clustal_path):
    '''
    Replaces the header line of a clustal alignment from a nonstandard
    variation like LocARNA's clustal header to the standard header
    found in constants.py.
    '''
    clustal_standard_header = 'CLUSTAL 2.0.10 multiple sequence alignment\n\n'

    clustal_lines = open(clustal_path).read().split('\n')
    assert clustal_lines[0][:len('CLUSTAL')] == 'CLUSTAL'
    clustal_lines[0] = clustal_standard_header
    open(clustal_path, 'w').write('\n'.join(clustal_lines))