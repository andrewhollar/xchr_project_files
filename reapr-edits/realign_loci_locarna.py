import os, glob, math, sys, subprocess, time, threading, random, argparse, traceback, multiprocessing, shutil
import utilities

def run_locarna(locarna, clustal_path, ungap_fasta_path, target_dir, target_file, max_diff, alifold_consensus_dp, guide_tree, verbose=False):
    """Run LocARNA.  Returns True if successful, or False otherwise"""

    if not os.path.isdir(target_dir): os.makedirs(target_dir)

    # Setup parameters
    acd_arg          = '--alifold-consensus-dp' if alifold_consensus_dp else ''
    guide_tree_arg   = '--treefile {0}'.format(guide_tree) if guide_tree else ''
    max_diff_arg     = '--max-diff {0}'.format(max_diff)
    ref_arg          = '--max-diff-aln={0}'.format(clustal_path)
    target_dir_arg   = '--tgtdir={0}'.format(target_dir)

    # Run LocARNA, (don't include --write-structure because it's unnecessary since intermediate directory will be deleted?)
    cmd = '%s %s %s %s %s --keep-sequence-order --verbose %s %s' \
        % (locarna, acd_arg, guide_tree_arg, ref_arg, max_diff_arg, target_dir_arg, ungap_fasta_path)

    start_time = time.time()
    # subprocess.Popen(cmd, shell=True, stdout=open('/dev/null','w'), stderr=subprocess.STDOUT).wait()
    subprocess.Popen(cmd, shell=True, stdout=subprocess.STDOUT, stderr=subprocess.STDOUT).wait()
    if verbose: print >>sys.stderr, cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds\n'

    # Copy final alignment 'results/result.aln' and delete
    # intermediate directory unless final alignment wasn't created,
    # probably due to impossible banding constraints
    result_path = os.path.join(target_dir, 'results', 'result.aln')
    if os.path.isfile(result_path):
        
        utilities.fix_clustal_header(result_path)  # Fix clustal header (remove non-standard LocARNA header)
        shutil.copyfile(result_path, target_file)   # Copy final alignment
        shutil.rmtree(target_dir)                   # Delete intermediate directory

        return True
    else:
        if verbose: print >>sys.stderr, 'Error: no alignment could be produced.'
        return False

def run_locarna_pool(x):
    try:
        return run_locarna(*x)
    except KeyboardInterrupt:
        pass
    except:
        print >>sys.stderr, traceback.print_exc()
