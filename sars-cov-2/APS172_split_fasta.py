#!/usr/bin/env python3
import os
import numpy as np

dir = "/Users/asherpreskasteinberg/Desktop/code/recombo/APS172_SC2_analysis/gapped"
file = os.path.join(dir, "sample_late-focus_aligned_gaps.fasta")

allseqs = open(file, "r")

names = []
positions = []
seqline = 0
for position, seq in enumerate(allseqs):
    if seq.startswith(">"):
        terms = seq.split(">")
        name = terms[len(terms)-1]
        name = name.replace("/", "_")
        names.append(str.rstrip(name))
        positions.append(position)

with open(file, "r") as master:
    master_full = master.readlines()

for name in names:
    fastafile = os.path.join(dir,"late-focus", name+".fasta")
    with open(fastafile, "w+") as fasta:
        fasta.write(">"+name+"\n")
        start = positions[0]
        end = positions[1]
        linenums = np.arange(start+1, end)
        for i in linenums:
            seqline = master_full[i]
            seqline = seqline.replace("a", "A")
            seqline = seqline.replace("t", "T")
            seqline = seqline.replace("c", "C")
            seqline = seqline.replace("g", "G")
            fasta.write(seqline)

genomelist = os.path.join(dir, "late-focus_list")
with open(genomelist, "w+") as genomes:
    for name in names:
        genomes.write(name+"\n")
