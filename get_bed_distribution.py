import os, sys

contig_dict = {}

with open(sys.argv[1], "r") as bed_in:
    for line in bed_in.readlines():
        line_tokens = line.split()
        # entry_length = int(line_tokens[2]) - int(line_tokens[1])
        entry_contig = line_tokens[0]
        if entry_contig in contig_dict:
            contig_dict[entry_contig] +=1
        else:
            contig_dict[entry_contig] = 1


for key in sorted(contig_dict.keys()):
    print(key, contig_dict[key])
