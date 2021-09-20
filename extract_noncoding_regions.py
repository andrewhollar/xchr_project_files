import sys


# Method to pad an integer with zeros on the left, this returns a string of length num_positions.
def pad_int(input_int, num_positions):
    str_rep = str(input_int)
    while len(str_rep) < num_positions:
        str_rep = "0" + str_rep
    return str_rep

gencode_gtf = "/home/ahollar/bed_files/gencode.v38.annotation.xchr.havana.gtf"
out_file = sys.argv[1]

previous_gene_end = 0

noncoding_regions = []
nc_id = 0

for line in open(gencode_gtf, "r").readlines():
    line_tokens = line.split()
    if line_tokens[2] == "gene":

        if int(line_tokens[3]) > int(previous_gene_end):
            nc_entry = ["chrX", str(previous_gene_end), str(line_tokens[3]), "nc_{}".format(pad_int(nc_id,8))]
            noncoding_regions.append(nc_entry)
            nc_id += 1
            
        previous_gene_end = line_tokens[4]



with open(out_file, "w") as nc_bed:
    for region in noncoding_regions:
        nc_bed.write(("\t".join(region) + "\n"))
