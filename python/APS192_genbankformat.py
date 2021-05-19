#!/usr/bin/env python3
import os
import numpy as np

dir = "/Users/asherpreskasteinberg/Desktop/code/recombo/APS192_HP_analysis"
file = os.path.join(dir, "doi_10.5061_dryad.8qp4n__v1/BIGSdb_gene-by-gene_alignment.xmfa")

allseqs = open(file, "r")

gene_names = []
strain_names = []
genome_positions = []
positions = []
seqline = 0
for position, seq in enumerate(allseqs):
    if seq.startswith(">"):
        ##get gene name
        terms = seq.split(">")
        header = terms[len(terms)-1]
        terms = header.split(" ")
        gene_name = terms[len(terms)-1]
        ##get strain name and position
        strain_pos = terms[0].split(":")
        strain = strain_pos[0]
        genome_position = strain_pos[1]
        gene_names.append(str.rstrip(gene_name))
        strain_names.append(str.rstrip(strain))
        genome_positions.append(str.rstrip(genome_position))
        positions.append(position)

with open(file, "r") as master:
    master_full = master.readlines()

outfile = os.path.join(dir, "APS192_HP_MASTER_MSA")
with open(outfile, "w+") as msa:
    for i in np.arange(0, len(positions)):
        msa.write(">"+gene_names[i]+" "+genome_positions[i]+" "+strain_names[i]+"\n")
        start = positions[i]
        if i == len(positions) - 1:
            end = len(master_full)
        else:
            end = positions[i+1]
        linenums = np.arange(start+1, end)
        for j in linenums:
            seqline = master_full[j]
            msa.write(seqline)

genome_names = set(strain_names)
genomelist = os.path.join(dir, "strain_list")
with open(genomelist, "w+") as genomes:
    for name in genome_names:
        genomes.write(name+"\n")
print(len(genome_names))