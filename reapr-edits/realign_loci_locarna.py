import os, glob, math, sys, subprocess, time, threading, random, argparse, traceback, multiprocessing, shutil
import utilities

#def run_locarna(locarna, clustal_path, ungap_fasta_path, target_dir, target_file, max_diff, alifold_consensus_dp, guide_tree, verbose=False):
def run_locarna(locarna, ungap_fasta_path, target_dir, target_file, max_diff, alifold_consensus_dp, guide_tree, verbose=False):

    """Run LocARNA.  Returns True if successful, or False otherwise"""
    
    run_locarna_flag = True
    
    log = []

    if not os.path.isdir(target_dir): os.makedirs(target_dir)

    # Setup parameters
    acd_arg          = '--alifold-consensus-dp' if alifold_consensus_dp else ''
    
    # -------------------------------------------------------------------------------
    # EDIT: The treefile argument now uses an '=' : --treefile=file
    guide_tree_arg   = '--treefile={0}'.format(guide_tree) if guide_tree else ''
    # -------------------------------------------------------------------------------
    # -------------------------------------------------------------------------------
    # EDIT: The max-diff argument now uses an '=' : --max-diff=difference
    max_diff_arg     = '--max-diff={0}'.format(max_diff)
    # ref_arg          = '--max-diff-aln={0}'.format(clustal_path)
    target_dir_arg   = '--tgtdir={0}'.format(target_dir)

    # Run LocARNA, (don't include --write-structure because it's unnecessary since intermediate directory will be deleted?)
    
    # -------------------------------------------------------------------------------
    # Edit: Updated the input parameters to match the Thiel pipeline.
    cmd = '%s --probabilistic --consistency-transformation --mea-beta=400 --iterations=2 %s %s %s --keep-sequence-order %s %s' \
        % (locarna, acd_arg, guide_tree_arg, max_diff_arg, target_dir_arg, ungap_fasta_path)
      #  % (locarna, acd_arg, guide_tree_arg, ref_arg, max_diff_arg, target_dir_arg, ungap_fasta_path)
    # -------------------------------------------------------------------------------

    locarna_results_dir = os.path.join(target_dir, "results", "")
    result_path = os.path.join(locarna_results_dir, 'result.aln')
    filtered_path = os.path.join(target_dir, 'improved_boundaries.aln')
    
    if os.path.isdir(locarna_results_dir):
        # locarna_result_aln = os.path.join(locarna_results_dir, "result.aln")
        locarna_reliabilities = os.path.join(locarna_results_dir, "result.bmreliability")
        
        if os.path.isfile(result_path) and os.stat(result_path).st_size != 0 and os.path.isfile(locarna_reliabilities) and os.stat(locarna_reliabilities).st_size != 0:
            if os.path.isfile(filtered_path) and os.stat(filtered_path).st_size != 0:
                # Locarna has already been run successfully.
                run_locarna_flag = False
            else:
                reliability_fit_once_out_path, reliability_log = run_reliability_profile_on_locarna_output(target_dir)
                log.extend(reliability_log)
                write_improved_boundaries_alignment(reliability_fit_once_out_path, result_path, filtered_path)
                run_locarna_flag = False
            
    if run_locarna_flag:
        start_time = time.time()
        # subprocess.Popen(cmd, shell=True, stdout=open('/dev/null','w'), stderr=subprocess.STDOUT).wait()
        # -------------------------------------------------------------------------------
        locarna_out_path = os.path.join(target_dir, "locarna.output")
        locarna_output = open(locarna_out_path, 'w', 1000000)    
        subprocess.Popen(cmd, shell=True, stdout=locarna_output, stderr=subprocess.STDOUT).wait()
        locarna_output.close()
        # -------------------------------------------------------------------------------
        
        #if verbose: print cmd + '\nRunning time: ' + str(time.time() - start_time) + ' seconds\n'
        log.append(cmd)
        time_str = 'Running time: ' + str(time.time() - start_time) + ' seconds'
        log.append(time_str)

        # Copy final alignment 'results/result.aln' and delete
        # intermediate directory unless final alignment wasn't created,
        # # probably due to impossible banding constraints
        # result_path = os.path.join(target_dir, 'results', 'result.aln')
        # filtered_path = os.path.join(target_dir, 'improved_boundaries.aln')
        
        # -------------------------------------------------------------------------------
        # EDIT: add run of reliability-profile.pl script to extract the high-confidence portion of the alignment
        reliability_fit_once_out_path, reliability_log = run_reliability_profile_on_locarna_output(target_dir)
        log.extend(reliability_log)
        write_improved_boundaries_alignment(reliability_fit_once_out_path, result_path, filtered_path)
    # -------------------------------------------------------------------------------

    if os.path.isfile(filtered_path) and os.stat(filtered_path).st_size != 0:        
        # utilities.fix_clustal_header(result_path)  # Fix clustal header (remove non-standard LocARNA header)
        # shutil.copyfile(result_path, target_file)   # Copy final alignment
        # shutil.rmtree(target_dir)                   # Delete intermediate directory

        return True, log
    else:
        # if verbose: print 'Error: no alignment could be produced.'
        log.append('Error: no alignment could be produced.')
        return False, log

def write_improved_boundaries_alignment(reliability_profile_output, result_alignment_path, filtered_result_path):
    reliability_profile_fit_once_output = open(reliability_profile_output).read().split('\n')
    if "Cannot read" not in reliability_profile_fit_once_output[0]:
        
        fit_line = ""
        for line in reliability_profile_fit_once_output:
            if "FIT" in line:
                fit_line = line            
            
        try:
            assert "FIT" in fit_line
        except AssertionError:
            print reliability_profile_output
            print result_alignment_path
            print filtered_result_path
            print fit_line
            print open(reliability_profile_output).read().split('\n')
            # raise AssertionError("End")
        
        
        fit_line_tokens = fit_line.split()
        assert fit_line_tokens[0] == "FIT"
        fit_line_tokens = fit_line_tokens[1:]
        
        fit_line_indices = [int(x) for x in fit_line_tokens]
        boundaries_start = min(fit_line_indices)
        boundaries_end = max(fit_line_indices)

        # boundaries_start = int(fit_line[1])
        # boundaries_end = int(fit_line[2])
        
        curr_align = open(result_alignment_path).read()
        assert curr_align != '', 'Alignment is empty'
        
        clustal_list = [x.split() for x in curr_align.split('\n') if x!='' and x[:7]!='CLUSTAL']
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

        # seq_length = int(seq_lengths[0])

        # for species, sequence in zip(header_list, seq_list):
        
        clustal_string = utilities.generate_clustal_boundaries(header_list, seq_list, boundaries_start, boundaries_end)
        open(filtered_result_path, 'w').write(clustal_string)

def run_reliability_profile_on_locarna_output(target_dir, fit_once_on = True):
    from commands import RELIABILITY_PROFILE
    
    log = []
    
    fit_once_on = '--fit-once-on' if fit_once_on else ''
    cmd = 'perl %s --dont-plot %s %s' % (RELIABILITY_PROFILE, fit_once_on, target_dir)

    start_time = time.time()
    # subprocess.Popen(cmd, shell=True, stdout=open('/dev/null','w'), stderr=subprocess.STDOUT).wait()
    # -------------------------------------------------------------------------------
    reliability_fit_once_out_path = os.path.join(target_dir, "reliability_fit_once.out")
    reliability_fit_once_output = open(reliability_fit_once_out_path, 'w', 1000000)    
    subprocess.Popen(cmd, shell=True, stdout=reliability_fit_once_output, stderr=subprocess.STDOUT).wait()
    reliability_fit_once_output.close()
    
    log.append(cmd)
    time_str = 'Running time: ' + str(time.time() - start_time) + ' seconds'
    log.append(time_str)
    
    # fit_once_on = ''
    # cmd = 'perl %s --dont-plot %s %s' % (RELIABILITY_PROFILE, fit_once_on, target_dir)
    # start_time = time.time()
    # reliability_out_path = os.path.join(target_dir, "reliability.out")
    # reliability_output = open(reliability_out_path, 'w', 1000000)    
    # subprocess.Popen(cmd, shell=True, stdout=reliability_output, stderr=subprocess.STDOUT).wait()
    # reliability_output.close()
    
    # log.append(cmd)
    # time_str = 'Running time: ' + str(time.time() - start_time) + ' seconds'
    # log.append(time_str)

    return reliability_fit_once_out_path, log
    
    # extracted_seq = open(extracted_output).read().split('\n')[1].strip()


def run_locarna_pool(x):
    try:
        return run_locarna(*x)
    except KeyboardInterrupt:
        pass
    except:
        print traceback.print_exc()
