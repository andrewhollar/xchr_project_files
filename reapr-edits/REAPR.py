import os, sys, multiprocessing, argparse, shutil, traceback, subprocess

import run_RNAz_screen
import tabulate_rnaz_results
import extract_loci
import realign_loci_locarna
import combine_tables
import commands
import utilities
import random


def write_REAPR_output(reapr_out):
    # PRINT THE OUTPUT TO OUT.log
    line = reapr_out.pop(0)
    try:
        while line:
            print line + '\n'
            line = reapr_out.pop(0)
    except IndexError:
        return

def main():
    REAPR_OUT_LINES = []
    
    parser = argparse.ArgumentParser(description='Realign a WGA and predict structural ncRNAs.')
    parser.add_argument('-a', '--alignments', required=True, help='Space-separated list of WGA alignment block files.')
    parser.add_argument('-s', '--species', required=True, help='Space-separated list of species in WGA.  Species names must be the same as those listed in the alignment block files.')
    parser.add_argument('-g', '--guide-tree', required=True, help='Species guide tree (in Newick format, without branch lengths) for progressive alignment by LocARNA.')
    parser.add_argument('-o', '--output-folder', default=os.getcwd(), help='Directory to write output files.  (Default: present working directory)')
    parser.add_argument('-d', '--delta', type=int, default=[20], nargs='+', help='Space-separated list of realignment deviations. (Default: 20)')
    parser.add_argument('-t', '--threshold', type=float, default=utilities.STABILITY_THRESHOLD, help='Stability filter threshold.  Filter out windows whose mean MFE z-score is above this threshold (Default: -1)')
    parser.add_argument('-p', '--processes', type=int, default=utilities.PROCESSES, help='Number of cores to use for multiprocessing (Default: 1)')
    parser.add_argument('-r', '--ram-disk', help='Location of RAM Disk to write temporary files.  Minimizes random access on disk storage.  This is highly recommended as REAPR will write many small files. (Note: this is typically /dev/shm in Ubuntu, and other Linux systems)')
    parser.add_argument('--alistat', action='store_true', help='Compute sequence identities of alignments using alistat')
    parser.add_argument('--compalignp', action='store_true', help='Compute change between original alignment and realignment using compalignp')
    args = parser.parse_args()

    # Setup temporary directory
    if args.ram_disk is None: tmp_dir = None
    else:
        utilities.directory_exists(args.ram_disk)
        tmp_dir = os.path.join(args.ram_disk, 'REAPR.'+utilities.get_random_string(8))
        os.makedirs(tmp_dir)
    
    NUM_PROCESSES = args.processes
    if NUM_PROCESSES == -1:
        NUM_PROCESSES = multiprocessing.cpu_count() - 2
    
    try:
        # Check if the MAF alignment block-index file exists.
        utilities.file_exists(args.alignments)
        
        # Load the MAF alignment block information into a list. Each entry in this list is a list containing two values:
        #       1: Name of the alignment block
        #       2: Filepath of the alignment block
        blocks = [x.split('\t') for x in open(args.alignments).read().split('\n') if x!='']
        block_dict = dict(blocks)
        block_names, block_paths = zip(*blocks)
        
        # Check if all of the alignment block filepaths exist.
        for a in block_paths: utilities.file_exists(a)

        # Check if the file containing the list of species exists
        utilities.file_exists(args.species)
        
        # Load the names of the species included in the MSA into a list.
        species = sorted([x for x in open(args.species).read().split('\n') if x!=''])
        
        ### Run RNAz screen on WGA ###
        
        # Locate RNAz and rnazWindow.pl
        # -------------------------------------------------------------------------------
        # EDIT: Changed the value of 'no_reference' to False, because we are using Human as a reference.    
        no_reference, structural, verbose, both_strands = False, False, True, True
        # -------------------------------------------------------------------------------
        alignment_format='MAF'
        
        # -------------------------------------------------------------------------------
        # EDIT: Changed the name of the output folder to match previous changes. All output
        #       will go into the 'alignments/' directory within the args.output_folder.
        out_dir = os.path.join(args.output_folder, 'alignments')
        # -------------------------------------------------------------------------------
        # -------------------------------------------------------------------------------
        # EDIT: Added a call to make the output directory to avoid downstream IOErrors
        if not os.path.isdir(out_dir): 
            os.makedirs(out_dir)
            
        # EDIT: Added directory to hold the reverse-complemented contigs to avoid doing redundant work.
        if not os.path.isdir(utilities.REV_COMP_CONTIG_DIR):
            os.makedirs(utilities.REV_COMP_CONTIG_DIR)
        # -------------------------------------------------------------------------------
        # -------------------------------------------------------------------------------
        # EDIT: Added the rnazSelectSeqs.pl as an argument.
        rnaz_1_args = [(alignment, no_reference, both_strands, utilities.WINDOW_SIZE, utilities.WINDOW_SLIDE, structural, commands.RNAz, commands.rnazWindow, commands.rnazSelectSeqs, out_dir, tmp_dir, alignment_format, 1, verbose) for alignment in block_paths] 
        # -------------------------------------------------------------------------------

        REAPR_OUT_LINES.append('Start: RNAz screen on WGA %s' % (utilities.get_time()))
        RNAZ_OUT_LINES = []
        
        pool = multiprocessing.Pool(processes=NUM_PROCESSES)        
        r = pool.map_async(run_RNAz_screen.eval_alignment_multiprocessing, rnaz_1_args, callback=RNAZ_OUT_LINES.extend) #.get(99999999)
        r.wait()
        for rnaz_entry in RNAZ_OUT_LINES:
            REAPR_OUT_LINES.extend(rnaz_entry)        
        REAPR_OUT_LINES.append('End: RNAz screen on WGA %s' % (utilities.get_time()))

        write_REAPR_output(REAPR_OUT_LINES)
        # print REAPR_OUT_LINES

        ### Compile table of RNAz screen results ###
        rnaz_paths = [a + '.rnaz' for a in block_paths]              # RNAz output
        log_paths = [a + '.windows.log' for a in block_paths]        # rnazWindow verbose logs
        index_paths = [a + '.windows.indices' for a in block_paths]  # window to slice indices map
        initial_table = os.path.join(args.output_folder, 'original_wga.tab')
        alternate_strands, merge = True, True
        tabulate_rnaz_results.write_table(initial_table, rnaz_paths, block_names, log_paths, index_paths, alternate_strands, merge, args.threshold, species)

        ### Extract stable loci ###
        loci_dir = os.path.join(args.output_folder, 'loci')
        
        # -------------------------------------------------------------------------------
        # EDIT: Changed the second to last argument to True, this indicates that the chromosome
        #       names need to be removed from the MAF species text.
        #       This list now contains tuples of length 2, (not 3) as I have removed the clustal.
        loci_alignment_list = extract_loci.extract_loci(block_dict, initial_table, args.threshold, loci_dir, species, utilities.WINDOW_SIZE, utilities.WINDOW_SLIDE, True, stdout=False)
        # -------------------------------------------------------------------------------

        realign_tables = [os.path.join(args.output_folder, 'locarna.d_%s.tab' % str(d)) for d in args.delta]

        for delta, realign_table in zip(args.delta, realign_tables):

            ### Realign loci ###

            # Setup alignment file paths 
            # old: locus_names, ref_clustals, ungap_fastas = zip(*loci_alignment_list)
            locus_names, ungap_fastas = zip(*loci_alignment_list)
            suffix = 'locarna{0}.{1}'.format('.g' if args.guide_tree else '', delta)
            target_dirs =  [x + '.%s.d' % suffix for x in ungap_fastas] #ref_clustals]
            # -------------------------------------------------------------------------------
            # EDIT: Changed the path of the target clustal file to be the boundaries predicted by the reliability-profile.pl script
            target_files = [os.path.join(x + '.%s.d'   % suffix, 'improved_boundaries.aln') for x in ungap_fastas] #ref_clustals]
            # -------------------------------------------------------------------------------


            # Realign loci
            acd, verbose = True, True
            locarna_target_args = [(commands.mlocarna, a, b, c, delta, acd, args.guide_tree, verbose) for a,b,c in zip(ungap_fastas, target_dirs, target_files)]

            REAPR_OUT_LINES.append('Start: LocARNA realignment, Delta=%s, %s' % (str(delta), utilities.get_time()))
            pool = multiprocessing.Pool(processes=NUM_PROCESSES)     
            
            LOCARNA_OUT_LINES = []     
            locarna_success = []   
            r = pool.map_async(realign_loci_locarna.run_locarna_pool, locarna_target_args, callback=LOCARNA_OUT_LINES.extend)
            r.wait()
                        
            for locarna_entry in LOCARNA_OUT_LINES:
                REAPR_OUT_LINES.extend(locarna_entry[1])
                locarna_success.append(locarna_entry[0])
                            
            target_files = [x for x,y in zip(target_files, locarna_success) if y]
            REAPR_OUT_LINES.append('End: LocARNA realignment, Delta=%s, %s' % (str(delta), utilities.get_time()))

            write_REAPR_output(REAPR_OUT_LINES)

            # # PRINT THE OUTPUT TO OUT.log
            # line = REAPR_OUT_LINES.pop(0)
            # try:
            #     while line:
            #         print line + '\n'
            #         line = REAPR_OUT_LINES.pop(0)
            # except IndexError:
            #     pass


            ### Run RNAz screen on realigned loci ###
            # -------------------------------------------------------------------------------
            # EDIT: Changed the value of 'no_reference' to False, because we are using Human as a reference. 
            #       Changed the value of 'both_strands' to True. 
            no_reference, structural, verbose, both_strands = False, True, True, True
            # -------------------------------------------------------------------------------
            alignment_format='CLUSTAL'
            rnaz_2_args = [(alignment, no_reference, both_strands, utilities.WINDOW_SIZE, utilities.WINDOW_SLIDE, structural, commands.RNAz, commands.rnazWindow, commands.rnazSelectSeqs, location, None, alignment_format, 2, verbose) for alignment, location in zip(target_files, target_dirs)]
            REAPR_OUT_LINES.append('Start: RNAz screen on realigned loci, Delta=%s, %s' % (str(delta), utilities.get_time()))
            
            RNAZ_OUT_LINES = []
            pool = multiprocessing.Pool(processes=NUM_PROCESSES)     
            r = pool.map_async(run_RNAz_screen.eval_alignment_multiprocessing, rnaz_2_args, callback = RNAZ_OUT_LINES.extend)
            r.wait()
            
            for rnaz_entry in RNAZ_OUT_LINES:
                REAPR_OUT_LINES.extend(rnaz_entry)
            
            REAPR_OUT_LINES.append('End: RNAz screen on realigned loci, Delta=%s, %s' % (str(delta), utilities.get_time()))
            
            # # PRINT THE OUTPUT TO OUT.log
            # line = REAPR_OUT_LINES.pop(0)
            # try:
            #     while line:
            #         print line + '\n'
            #         line = REAPR_OUT_LINES.pop(0)
            # except IndexError:
            #     pass
            write_REAPR_output(REAPR_OUT_LINES)

            
            ### Compile tables of RNAz results ###
            rnaz_paths = [a + '.rnaz' for a in target_files]              # RNAz output
            
            # -------------------------------------------------------------------------------
            # EDIT: There are no windows in the second pass of RNAz
            # log_paths = [a + '.windows.log' for a in target_files]        # rnazWindow verbose logs
            log_paths = [None for a in target_files]
            #index_paths = [a + '.windows.indices' for a in target_files]  # window to slice indices map
            index_paths = [None for a in target_files]
            # -------------------------------------------------------------------------------
            alternate_strands, merge = True, False
            block_names = locus_names
            tabulate_rnaz_results.write_table(realign_table, rnaz_paths, block_names, log_paths, index_paths, alternate_strands, merge, args.threshold, species)

        # Combine tables
        combine_tables.combine_tables(initial_table, realign_tables, args.delta, loci_dir, True, species, args.alistat, args.compalignp, os.path.join(args.output_folder, 'summary.tab'))

    except KeyboardInterrupt:
        print 'Outside except, terminating pool'
        pool.terminate()

    except:
        print traceback.print_exc()
    finally:
        # # Kill descendant processes, i.e. those with same group id
        # group_processes = subprocess.Popen('pgrep -g %s' % (parent_pid), shell=True, stdout=subprocess.PIPE).communicate()[0].split()
        # group_processes.remove(str(parent_pid))
        # subprocess.Popen('kill -9 %s' % (' '.join(group_processes)), shell=True).wait()

        # Clear temporary directory in RAM Disk
        if tmp_dir and os.path.isdir(tmp_dir): shutil.rmtree(tmp_dir)


if __name__=='__main__':
    main()
