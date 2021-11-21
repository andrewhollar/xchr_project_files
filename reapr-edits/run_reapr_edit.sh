#!/bin/bash
OUTDIR=/home/ahollar/reapr_x_full

#This is where the MAF preprocessing is done (Python 3)
# source ~/mafPreprocessingEnv/bin/activate
# python ~/xchr_project_files/reapr-edits/process_maf_file.py ~/alignments/7way_final.maf $OUTDIR
# deactivate

#This is where REAPR is run (Python 2)
source ~/reaprEnv/bin/activate
python ~/xchr_project_files/reapr-edits/REAPR.py -p -1 -a $OUTDIR/alignment_blocks.txt -s ~/7way.species -g ~/7way.newick -o $OUTDIR -t 0.5 > ~/reapr_x_full2.log 
deactivate