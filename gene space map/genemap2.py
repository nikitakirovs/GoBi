import numpy as np
import csv
import glob

unique_genome = []
gaps = [
    [6, 11, 57, 0],
    [3, 49, 54, 0],
    [3, 14, 26, 0],
    [5, 19, 0],
    [1, 0, 0],
    [1, 0, 0],
    [1, 0, 0]
]

matrix_colors = [
    [0, 800, 800, 3400, 3400, 3700],
    [2900, 3700, 300, 2900, 0, 300],
    [300, 1700, 1700, 3000, 0, 300],
    [400, 1200, 0, 400, 1200, 3000],
    [200, 3000, 0, 0, 0, 200],
    [250, 3000, 0, 0, 0, 250],
    [1, 1, 1, 1, 1, 3000],
]

# get all the regions of reference(Apis melliefera)
matrix = np.array([1, 2, 3])  # genome, start, end
with open('regions/Nasonia_vitripennis.txt', mode ='r')as file:
  csvFile = csv.reader(file)
  for lines in csvFile:
      line = lines[0].split('landmark:')[1]
      line = line.split('_')
      genome = "Apis_mellifera"
      if genome not in unique_genome:
          unique_genome.append(genome)
      start = line[4].split(':')[1]
      end = line[5].split(':')[1]
      new_row = np.array([genome, start, end])
      matrix = np.vstack([matrix, new_row])

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

def filter_rows_with_same_first_element(element, matrix):
    filtered_rows = []
    for row in matrix:
        if row[0] == element:
            filtered_rows.append(row)
    return filtered_rows

def get_genes_list(element, matrix, number):
    k = 0
    same_first_element_rows = filter_rows_with_same_first_element(element, matrix)
    list_of_genes = ()
    vstart = 100
    vend = 140
    for i in range(len(same_first_element_rows)):
        start = int(same_first_element_rows[i][1])
        end = int(same_first_element_rows[i][2])
        if(start < end):
            strand = 1
        else:
            strand = -1
        new_gene = (vstart, vend, strand)
        if(gaps[number][k] - 1 == i):
            gap = 100
            k = k + 1
        else:
            gap = 15

        list_of_genes = list_of_genes + (new_gene,)
        vstart = vstart + 40 + gap
        vend = vend + 40 + gap

    return list_of_genes

genome_list = [{1:100}]

import matplotlib.colors as colors
colors_list = list(colors._colors_full_map.values())

number = 0
myorder = [0, 5, 1, 3, 4, 6, 2]
unique_genome = [unique_genome[i] for i in myorder]
for genes in unique_genome:
    new_line = {"name": genes, "size": 3700, "cds_list": get_genes_list(genes, matrix, number)}
    genome_list.append(new_line)
    number = number + 1

from pygenomeviz import GenomeViz

gv = GenomeViz(fig_track_height = 0.5)
colors = ("skyblue", "tomato", "green", "grey")

for i in range(1, len(genome_list)):
    name, size, cds_list = genome_list[i].get("name"), genome_list[i].get("size"), genome_list[i].get("cds_list")
    track = gv.add_feature_track(name, size)
    j = 0
    for idx, cds in enumerate(cds_list, 1):
        size = len(cds_list)
        start, end, strand = cds
        point1 = matrix_colors[i - 1][0]
        point2 = matrix_colors[i - 1][1]
        point3 = matrix_colors[i - 1][2]
        point4 = matrix_colors[i - 1][3]
        point5 = matrix_colors[i - 1][4]
        point6 = matrix_colors[i - 1][5]
        color = 0
        if start > point1 and end < point2:
            color = 0
        if start > point3 and end < point4:
            color = 1
        if start > point5 and end < point6:
            color = 2

        track.add_feature(start, end, strand, facecolor=colors[color])
        j = j + 1

gv.savefig("genemap2.png")