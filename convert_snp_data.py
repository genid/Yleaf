
hg19_data = []
hg38_data = []
with open("snp_table.csv") as f:
    for line in f:
        try:
            name, branch, hg19_location, hg38_location, ancestor, derived = line.strip().split(",")
        except ValueError:
            continue
        hg19_data.append(["chry", name, branch, hg19_location, f"{ancestor}->{derived}", ancestor, derived])
        hg38_data.append(["chry", name, branch, hg38_location, f"{ancestor}->{derived}", ancestor, derived])

with open("../Position_files/WGS_hg19_new.txt", "w") as f:
    for line in hg19_data:
        f.write('\t'.join(line) + "\n")

with open("../Position_files/WGS_hg38_new.txt", "w") as f:
    for line in hg38_data:
        f.write('\t'.join(line) + "\n")
