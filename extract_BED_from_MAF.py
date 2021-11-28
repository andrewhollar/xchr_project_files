import os
import sys
import shutil
import argparse

from Bio import AlignIO

def run(fn, ref):
    maf = AlignIO.parse(fn, "maf")
    for msa in maf:
        for sr in msa:
            if sr.id.startswith(ref):
                contig = sr.id.split('.')[1]
                start_pos = sr.annotations['start']
                end_pos = int(sr.annotations['start']) + int(sr.annotations['size'])
                name = sr.id
                score = "0"
                strand = ""
                
                if sr.annotations['strand'] == -1:
                    strand = "-"
                else:
                    strand = "+"
                
                bed_line = "\t".join([contig, str(start_pos), str(end_pos), name, score, strand])
                print(bed_line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", '--maf-file',
                        type=str, required=True,
                        help="input MAF file")
    parser.add_argument("-r", '--reference-species',
                        type=str,
                        help="reference species, e.g., mm10, hg19")
    ARGS = parser.parse_args()
    MAF_FILEPATH = ARGS.maf_file
    REFERENCE_SPECIES = ARGS.reference_species
    try:
        run(MAF_FILEPATH, REFERENCE_SPECIES)
    except FileNotFoundError:
        print("File {} not found".format(MAF_FILEPATH))
        exit(1)
    exit(0)