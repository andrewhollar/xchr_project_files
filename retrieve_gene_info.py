import argparse
import os
from pyensembl import EnsemblRelease

UCSC_BED_DIR = "/home/ahollar/bed_files/ucsc_hg38_gencode"
UCSC_WHOLE_GENES = os.path.join(UCSC_BED_DIR, "ucsc_chrX_wholeGenes.bed")

MAX_SAMPLES = 5

def main():
    # release 77 uses human reference genome GRCh38
    data = EnsemblRelease(104)

    # gene_name = data.gene_name_of_transcript_id("ENST00000431238") 
    # gene_info = data.genes_by_name(gene_name)
    # print(gene_name)
    # print(gene_info)
    
    entries = 0
    
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

        # print(transcripts)
        
        if 'pseudogene' in transcript.biotype:
            print(transcript)
            entries +=1
        # print(transcript_id)
        
        
        # gene_name = data.gene_name_of_transcript_id(transcript_id.split('.')[0])
        
        # ids_for_gene = data.transcript_ids_of_gene_name(gene_name)
        # print(ids_for_gene)
        # gene_info = data.genes_by_name(gene_name)
        # print(gene_name)
        # print(gene_info)


        
        if entries > MAX_SAMPLES:
            break
        



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