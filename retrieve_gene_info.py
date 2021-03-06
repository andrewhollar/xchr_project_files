import argparse
import os
import sys
from pyensembl import EnsemblRelease

UCSC_BED_DIR = "/home/ahollar/bed_files/ucsc_hg38_gencode"
UCSC_WHOLE_GENES = os.path.join(UCSC_BED_DIR, "ucsc_chrX_wholeGenes.bed")

MAX_SAMPLES = sys.maxsize

def main():
    # release 77 uses human reference genome GRCh38
    data = EnsemblRelease(104)

    entries = 0
    
    biotype_cts = {}
    
    for entry in open(UCSC_WHOLE_GENES, "r").readlines():
        entry_tokens = entry.split()
        start_pos = entry_tokens[1]
        end_pos = entry_tokens[2]
        transcript_id = entry_tokens[3]
        
        try:
            transcript = data.transcript_by_id(transcript_id.split('.')[0])
        except ValueError:
            print("ValueError: Transcript not found: {}".format(transcript_id.split('.')[0]))
            continue        
        
        if transcript.biotype in biotype_cts.keys():
            biotype_cts[transcript.biotype] += 1
        else:
            biotype_cts[transcript.biotype] = 1
        
        if entries > MAX_SAMPLES:
            break
        
    print(biotype_cts)



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
        #run(MAF_FILEPATH, REFERENCE_SPECIES)
    except FileNotFoundError:
        # print("File {} not found".format(MAF_FILEPATH))
        exit(1)
    exit(0)