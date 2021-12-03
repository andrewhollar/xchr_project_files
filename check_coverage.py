import os
import sys

UCSC_BED_DIR = "/home/ahollar/bed_files/ucsc_hg38_gencode"
UCSC_EXONS = os.path.join(UCSC_BED_DIR, "ucsc_chrX_exons.bed")
UCSC_INTRONS = os.path.join(UCSC_BED_DIR, "ucsc_chrX_introns.bed")
UCSC_INTERGENIC = os.path.join(UCSC_BED_DIR, "ucsc_chrX_intergenic.bed")

MAX_SAMPLES = sys.maxsize

def get_num_nts(bed_file):
    num_nts = 0
    for entry in open(bed_file, "r").readlines():
        entry_tokens = entry.split()
        start_pos = entry_tokens[1]
        end_pos = entry_tokens[2]
        #transcript_id = entry_tokens[3]
        assert (end_pos - start_pos) > 0
        num_nts += (end_pos - start_pos)

    return num_nts

def main():
    
    exon_cov = get_num_nts(UCSC_EXONS)
    intron_cov = get_num_nts(UCSC_INTRONS)
    intergenic_cov = get_num_nts(UCSC_INTERGENIC)
    
    total = (exon_cov + intron_cov + intergenic_cov)
    
    print(total)
    print("Exon: {}".format((exon_cov/total)))
    print("Intron: {}".format((intron_cov/total)))
    print("Intergenic: {}".format((intergenic_cov/total)))



if __name__ == '__main__':
    # parser = argparse.ArgumentParser()
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
    except FileNotFoundError:
        exit(1)
    exit(0)