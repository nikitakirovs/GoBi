import numpy as np
import csv
import glob

#length of the regions which come before the cds of our interests on the genome
pre_length = (5760000, 100000, 230000, 2100000, 4660000, 66160000, 2230000)

motifs = dict()
#read fasta file of Apis mellifera and create the first row of gene map based on it
with open('regions/Apis-Mellifera-cds-reduced.fasta') as f:
    lines = f.readlines()
    for i in range(0, len(lines)):
        s = lines[i].strip()
        if s[0] == '>':
            key = s[1:]
        else:
            motifs[key] = s

unique_genome = []
matrix = np.array([1, 2, 3])

files = glob.glob("regions/*.txt")
for file in files:
    with open(file, mode='r') as file:
        csvFile = csv.reader(file)
        for lines in csvFile:
            line = lines[0].split('n/a')[1]
            line = line.split()
            genome = file.name.split("\\")[1].split(".")[0]
            if genome not in unique_genome:
                unique_genome.append(genome)
            start = line[1]
            end = line[2]
            new_row = np.array([genome, start, end])
            matrix = np.vstack([matrix, new_row])



# for fasta files
key_list = list(motifs.keys())
val_list = list(motifs.values())
Apis_matrix = np.array([1, 2, 3, 4])

for i in range(len(motifs)):
    split = key_list[i].split("_")
    genome = "Apis mellifera"
    start = key_list[i].split("start:", 1)[1].split("_")[0]
    end = key_list[i].split("end:", 1)[1].split("_")[0]
    product = key_list[i].split("product=", 1)[1].split(";")[0]
    new_row = np.array([genome, start, end, product])
    Apis_matrix = np.vstack([Apis_matrix, new_row])

def filter_rows_with_same_first_element(element, matrix):
    filtered_rows = []
    for row in matrix:
        if row[0] == element:
            filtered_rows.append(row)
    return filtered_rows

def get_genes_list(element, matrix, number):
    same_first_element_rows = filter_rows_with_same_first_element(element, matrix)
    list_of_genes = ()
    for i in range(len(same_first_element_rows)):
        start = int(same_first_element_rows[i][1]) - pre_length[number]
        end = int(same_first_element_rows[i][2]) - pre_length[number]
        if(start < end):
            strand = 1
        else:
            tmp = end
            end = start
            start = tmp
            strand = -1
        new_gene = (start, end, strand)
        if new_gene not in list_of_genes:
            list_of_genes =  list_of_genes + (new_gene,)
    return list_of_genes

genome_list = [{1:100}]

new_line = {"name": "Apis mellifera", "size": 123000, "cds_list": get_genes_list("Apis mellifera", Apis_matrix, 6)}
genome_list.append(new_line)

number = 0
for genes in unique_genome:
    new_line = {"name": genes, "size": 123000, "cds_list": get_genes_list(genes, matrix, number)}
    genome_list.append(new_line)
    number = number + 1

from pygenomeviz import GenomeViz

gv = GenomeViz(fig_track_height = 0.5)
colors = ("skyblue", "tomato", "green")
matrix_colors = [
    [0, 20, 20, 102, 103, 120],
    [40, 60, 60, 90, 20, 30],
    [1, 1, 1, 1, 1, 120],
    [40, 60, 20, 35, 60, 120],
    [40, 120, 0, 0, 20, 40],
    [60, 120, 26, 60, 0, 25],
    [40, 120, 0, 0, 0, 40]
]


for i in range(1, len(genome_list)):
    name, size, cds_list = genome_list[i].get("name"), genome_list[i].get("size"), genome_list[i].get("cds_list")
    track = gv.add_feature_track(name, size)
    j = 0
    for idx, cds in enumerate(cds_list, 1):
        size = len(cds_list)
        start, end, strand = cds
        point1 = matrix_colors[i - 1][0] * 1000
        point2 = matrix_colors[i - 1][1] * 1000
        point3 = matrix_colors[i - 1][2] * 1000
        point4 = matrix_colors[i - 1][3] * 1000
        point5 = matrix_colors[i - 1][4] * 1000
        point6 = matrix_colors[i - 1][5] * 1000
        color = 0
        if start > point1 and end < point2:
            color = 0
        if start > point3 and end < point4:
            color = 1
        if start > point5 and end < point6:
            color = 2

        if (name == "Apis mellifera"):
            track.add_feature(1000, 1002, 1, label="yellow-like", linewidth=1, labelrotation=0, labelvpos="top", labelhpos="center", labelha="center")
            track.add_feature(21000, 21001, 1, label="Mrjps", linewidth=1, labelrotation=0, labelvpos="top",
                              labelhpos="center", labelha="center")
            track.add_feature(108000, 108001, 1, label="yellow-h", linewidth=1, labelrotation=0, labelvpos="top",
                              labelhpos="center", labelha="center")

        track.add_feature(start, end, strand, facecolor=colors[color])
        j = j + 1


gv.savefig("genemap.png")
