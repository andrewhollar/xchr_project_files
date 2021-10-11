import os, sys

length_dict = {}

with open(sys.argv[1], "r") as bed_in:
    for line in bed_in.readlines():
        line_tokens = line.split()
        entry_length = int(line_tokens[2]) - int(line_tokens[1])

        if entry_length in length_dict:
            length_dict[entry_length] +=1
        else:
            length_dict[entry_length] = 1


for key in sorted(length_dict.keys()):
    print(key, length_dict[key])
