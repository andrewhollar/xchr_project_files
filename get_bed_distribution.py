"""
This script reads a BED file and counts the number of entries on each of the sequences included.
The results are printed out to the console. 

This script requires one input parameter -- the BED file to analyze.

example:
python get_bed_distribution.py /path/to/BED/file

Output:
chr 1   83575
chr 2   77347
...     ... 

"""

import sys

contig_dict = {}

# Open the BED file
with open(sys.argv[1], "r") as bed_in:

    # Read through BED file line by line. Each line in a BED file contains the entry for one genomic region.
    for line in bed_in.readlines():
        line_tokens = line.split()

        # Parse the sequence of the genomic region. This is the first column in the BED format.
        entry_contig = line_tokens[0]
        if entry_contig in contig_dict:
            contig_dict[entry_contig] +=1
        else:
            contig_dict[entry_contig] = 1

# Print the counts 
for key in sorted(contig_dict.keys()):
    print(key, contig_dict[key])
