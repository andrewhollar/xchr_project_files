import os, sys

with open(sys.argv[1], "r") as chromInfo:
    for line in chromInfo.readlines():
        line_tokens = line.split("|")

        if len(line_tokens) > 1:
            print(line_tokens)
            break
