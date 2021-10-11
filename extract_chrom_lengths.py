import os, sys

out_file = os.path.join(os.path.dirname(sys.argv[1]), "hg38.genome")

out_lines = []

with open(sys.argv[1], "r") as chromInfo:
    for line in chromInfo.readlines():
        line_tokens = line.split("|")

        if len(line_tokens) > 1:
            # print(line_tokens)

            out_lines.append("{}\t{}\n".format(line_tokens[1], line_tokens[2]))


with open(out_file, "w") as chromInfoBed:
    for line in out_lines:
        chromInfoBed.write(line)