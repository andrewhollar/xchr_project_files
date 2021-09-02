from Bio import AlignIO
import argparse



from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
# from config import CONF
from pathlib import Path
import sys
import re

annotation_tags = ["src", "start", "size", "strand", "srcSize"]
SMALLEST_SEQUENCE = 6


'''
This script is used to preprocess the 241-way MAF alignment. This scri
'''




#Extract the species included from the newick-formatted tree
#   newick_tree_path - the path to a plain text file containing a newick-formatted phylogenetic tree
#
#   returns a list containing strings of each species included in the phylogenetic tree
def extract_species_from_newick(newick_tree_path):

    #Read the text version of the tree
    nt = open(newick_tree_path, "r")
    newick_tree_str = nt.read()
    nt.close()

    #Remove all numeric/punctuation (except for ':' and '_'), all instances of 'fullTreeAnc' 
    newick_tree_str = newick_tree_str.translate(str.maketrans('', '', '()0123456789.,'))
    newick_tree_str = re.sub(r'fullTreeAnc', "", newick_tree_str)

    #Split species by the remaining ":" characters
    newick_tree_list = newick_tree_str.split(":")

    #Keep only the entries that are non-empty
    target_species = [species for species in newick_tree_list if species != ""][:-1]



    csl = ""
    for s in target_species:
        if s.startswith("e-"):
            print(s)
            s = s[2:]
            print(s)

        if s == "point":
            print("skipping entry point")
            continue

        for c in s:
            if c != " ":
                csl += c
        csl += ","
    print(csl)
    
    return target_species
    

def process_maf_file(arguments):

    MAX_ITERATIONS = 100

    try:
        iteration = 0
        for alignment_block in AlignIO.parse(arguments.maf_path, "maf"):
            if iteration > MAX_ITERATIONS:
                break
            print('processing alignment block')
            for species in alignment_block:
                print(species.id, species.annotations["start"], species.annotations["end"])

            iteration += 1
    except Exception as e:
        print(e)


# Function that takes a multiple alignment 
def main():
    parser = argparse.ArgumentParser(description='Preprocess the 241-way mammalian whole genome alignment (WGA).')
    parser.add_argument("-p", "--mafpath", help="The file path to the MAF-formatted alignment file.")
    parser.add_argument("--targetspecies", help="The species used as the reference when the MAF was created (e.g. Homo_sapiens).")
    parser.add_argument("--targetsequence", help="The name of the sequence to be extracted from the WGA (e.g. chrX).")
    parser.add_argument("--start", help="The start position on the TARGET_SEQUENCE of the region to extract.", type=int)
    parser.add_argument("--end", help="The end position on the TARGET_SEQUENCE of the region to extract.", type=int)


    process_maf_file(parser)

    #process_maf_file(sys_argv[1])

    #Read newick-formatted tree to determine which entries from MAF alignment blocks to remove (i.e. entries that contain 'fullTreeAnc' in the species name)
    #target_species = extract_species_from_newick(sys.argv[1])

    # print("Reading and filtering .maf file based on preferences set in config.py.")
    # print("The target species are: " + str(CONF['TARGET_SPECIES']))
    # m = MAF(Path(sys.argv[1]))
    # print("Filtering complete.")
    # print("Writing new .maf file")
    # m.write(Path(sys.argv[2]))
    # print("Write complete.")




# Class that represents a single alignment block within a .maf (multiple alignment format) file.
# This enables the score of the alignment block to be carried over from the unfiltered .maf to the
# updated, filtered .maf file. 
class MAFBlock:
    def __init__(self, score, seq_records = []):
        self.alignment = seq_records
        self.score = score
        self.length_flag = False

    # Add SeqRecord object to list of alignments
    def add_seq_record(self, sr):
        self.alignment.append(sr)

    # Yield all of the SeqRecords within the alignment
    def get_seq_records(self):
        for a in self.alignment:
            yield a

    # Filter the MAFBlock to only include species included as target species. Remove common gaps within
    # filtered alignment.
    def filter(self, in_alignment):
        filtered_srs = []
        for sr in self.alignment:
            id = sr.id.split('.')[0]
            if (id in CONF['TARGET_SPECIES']) or (id == CONF['REF_SPECIES']):
                #if in_alignment[CONF['TARGET_SPECIES'].index(id)]:
                filtered_srs.append(sr)

        self.alignment = filtered_srs
        self.remove_common_whitespace()

    # Get all characters at position idx of the alignment block, used to check for common gaps.
    def get_column(self, idx):
        col = ""
        for sr in self.get_seq_records():
            col += sr.seq[idx]
        return col

    # Remove a column from the alignment block, this occurs when all characters are gap-symbols 
    # i.e. there is a common gap
    def remove_column(self, idx):
        assert self.get_column(idx).count("-") == len(self.get_column(idx))

        for sr in self.alignment:
            sr.seq = sr.seq[:idx] + sr.seq[idx + 1:]

    # Loop through each of the characters within the alignment block and check for common gaps. 
    def remove_common_whitespace(self):
        seq_len = len(self.alignment[0].seq)
        i = 0

        while (i < seq_len):
            col_aln = self.get_column(i)
            if (col_aln.count("-") == len(col_aln)):
                self.remove_column(i)
                seq_len -= 1
            else:
                i += 1


# Class that represents a single .maf file. Comprised of a header, footer, and alignment blocks. The
# header and footer are both lists of strings, each coming from a comment line in the .maf file. The 
# alignment blocks are of type MAFBlock, to facilitate the transfer of scores.
class MAF:
    def __init__(self, mafPath = ""):
        self.blocks = []
        self.header = []
        self.footer = []
        if mafPath:
            self.read(mafPath)      

    # Add a new alignment MAFBlock to the list of blocks. 
    def add_block(self, block: MAFBlock):
        self.blocks.append(block)

    # Add a new header line to the header field.
    def add_header_line(self, line):
        self.header.append(line)

    # Add a new footer line to the footer field.
    def add_footer_line(self, line):
        self.footer.append(line)

    # Read/Parse a .maf file at the specified mafPath. Filter blocks as they are read. 
    def read(self, mafPath):
        current_block = MAFBlock(0.0, [])
        first = True

        for line in open(mafPath, "r"):
            # Encountered a comment line, either in the header or footer
            if line.startswith("#"):
                if first:
                    self.add_header_line(line.rstrip())
                else:
                    self.add_footer_line(line.rstrip())

            # Encountered a sequence within an alignment block
            elif line.startswith('s'):

                #Remove enpty strings from parsed line
                tokens = list(filter(None, line.split(' ')))
                
                #Information for annotation
                annotation_items = {}

                for t in range(1, len(tokens) - 1, 1):
                    annotation_items[annotation_tags[t - 1]] = tokens[t]
                    
                # Create new SeqRecord
                sr = SeqRecord(Seq(tokens[6].rstrip()), 
                               id = annotation_items['src'], 
                               annotations = annotation_items)

                print(tokens[3])
                if int(tokens[3]) < SMALLEST_SEQUENCE:
                    current_block.length_flag = True

                current_block.add_seq_record(sr)

            # Encountered a new alignment block
            elif line.startswith('a'):
                # This is the first alignment block that has been encountered
                if first:
                    first = False
                else:
                    # Check the existing alignment block for target species before overwriting it.
                    # If the reference and at least 1 other species are included, the MAFBlock will be
                    # filtered and saved.
                    if (not current_block.length_flag):
                        self.check_for_target_species(current_block)

                # Create new block using the score found in .maf file
                current_block = MAFBlock(line.split('=')[1].rstrip(), [])
            
            elif line.startswith('i'):
                continue

            elif line.startswith('e'):
                continue

            # Encountered an empty line
            elif line.rstrip() == "":
                continue
            else: 
                print(line)

        self.check_for_target_species(current_block)


    # Check the current MAFBlock for the target species declared in config.py. This 
    def check_for_target_species(self, block: MAFBlock):
        # Boolean list to keep track of which of the target species are included in the alignment block
        in_alignment = [False] * len(CONF['TARGET_SPECIES'])   

        # Loop through each of the SeqRecords within the alignment block, updating Boolean list as needed
        for sr in block.get_seq_records():
            id = sr.id.split('.')[0]
            if id in CONF['TARGET_SPECIES']:
                in_alignment[CONF['TARGET_SPECIES'].index(id)] = True   

        # If any of the target species are included in the alignment block, write this to new  
        if any(in_alignment[:]): 
            block.filter(in_alignment)
            
            self.add_block(block)
        
        # None of the target species were included in the alignment block (besides the reference)
        else:
            return


#  ["src", "start", "size", "strand", "srcSize"]
    def write(self, outPath):
        with open(outPath, "w") as f:
            for line in self.header:
                f.write(line + "\n")

            f.write("\n")
            
            for filteredBlock in self.blocks:
                f.write("a score=" + filteredBlock.score + "\n")
                for alignedSpecies in filteredBlock.get_seq_records():
                    lineString = "{0} {1:<12} {2:<6} {3:<6} {4} {5} {6}\n"
                    f.write(lineString.format("s", alignedSpecies.id, alignedSpecies.annotations['start'], alignedSpecies.annotations['size'], alignedSpecies.annotations['strand'], alignedSpecies.annotations['srcSize'], alignedSpecies.seq[:]))

                f.write("\n")

            f.write("\n")

            for line in self.footer:
                f.write(line + "\n")

if __name__ == "__main__":
    main()