#!/bin/bash

#cat myfile.txt | grep -n "12:09" | head -1 | cut -d":" -f1

CHR_X_START="$(grep -n -m1 'Homo_sapiens.chr1' ~/alignments/sample.maf | cut -d ':' -f1)"
CHR_X_END="$(grep -n 'Homo_sapiens.chr1' ~/alignments/sample.maf | tail -1 | cut -d ':' -f1)"
echo $CHR_X_START
echo $CHR_X_END


sed -n '${CHR_X_START},${CHR_X_END} p' /home/ahollar/alignments/sample.maf > extracted_region.maf