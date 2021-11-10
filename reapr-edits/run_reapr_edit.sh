#!/bin/bash

OUTDIR=/home/ahollar/reapr_x

#This is where the MAF preprocessing is done (Python 3)
source ~/mafPreprocessingEnv/bin/activate
python ~/xchr_project_files/reapr-edit/process_maf_file.py ~/alignments/6way_final.maf $OUTDIR
deactivate

#This is where REAPR is run (Python 2)
source ~/reaprEnv/bin/activate

deactivate