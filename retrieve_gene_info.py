import argparse
import os
from pyensembl import EnsemblRelease


def main():
    # release 77 uses human reference genome GRCh38
    data = EnsemblRelease(104)

    gene_info = data.gene_by_id("ENST00000431238.7")
    print(gene_info)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("-m", '--maf-file',
    #                     type=str, required=True,
    #                     help="input MAF file")
    # parser.add_argument("-r", '--reference-species',
    #                     type=str,
    #                     help="reference species, e.g., mm10, hg19")
    # ARGS = parser.parse_args()
    # MAF_FILEPATH = ARGS.maf_file
    # REFERENCE_SPECIES = ARGS.reference_species
    try:
        main()
        #run(MAF_FILEPATH, REFERENCE_SPECIES)
    except FileNotFoundError:
        # print("File {} not found".format(MAF_FILEPATH))
        exit(1)
    exit(0)