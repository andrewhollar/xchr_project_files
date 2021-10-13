import os, sys

out_file = os.path.join(os.path.dirname(sys.argv[1]), "hg38.genome")

out_lines = []

with open(sys.argv[1], "r") as chromInfo:
    for line in chromInfo.readlines():
        if "size" in line:
            continue

        line_tokens = line.split("|")

        if len(line_tokens) > 1:
            out_lines.append("{}\t{}\n".format(line_tokens[1].strip(), line_tokens[2].strip()))


with open(out_file, "w") as chromInfoBed:
    for line in out_lines:
        chromInfoBed.write(line)