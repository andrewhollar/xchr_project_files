#!/bin/bash

cd ~/alignments

mkdir 7way_genomes
cd 7way_genomes

GCA='GCA'
GCF='GCF'

cat ~/xchr_project_files/7way_accession_numbers.txt | while read -r asc ; do
	esearch -db assembly -query $asc </dev/null \
	| esummary \
	| xtract -pattern DocumentSummary -element FtpPath_GenBank \
	| while read  -r url ; do 
        if [[ "$asc" == *"$GCF"* ]]; then
            fname=$(echo $url | grep -o 'GCF_.*' | sed 's/$/_genomic.fna.gz/') ;
        else
		    fname=$(echo $url | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;
        fi
		wget "$url/$fname" ;
	done ;
done