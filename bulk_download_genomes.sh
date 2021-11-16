#!/bin/bash

cd ~/alignments

mkdir 7way_genomes
cd 7way_genomes

cat ~/xchr_project_files/7way_accession_numbers.txt | while read -r asc ; do
	esearch -db assembly -query $asc </dev/null \
	| esummary \
	| xtract -pattern DocumentSummary -element FtpPath_GenBank \
	| while read  -r url ; do 
		fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
		wget "$url/$fname" ;
	done ;
done