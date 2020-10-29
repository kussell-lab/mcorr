import sys
import argparse
import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
import pandas as pd
from tqdm import tqdm
import os

"""
Goal:

This code is for splitting MSA files containing all sequences into smaller MSA files
with equal numbers of sequences to be used with mcorr-pair-sync. This allows for parallel
computation of pairwise distances between sequences

"""

"""
this is pretty sweet for building help and parsing commandline arguments can be called
via python3 splitMSA.py -h
check out here https://docs.python.org/3/howto/argparse.html#id1
"""
parser = argparse.ArgumentParser(description='Splits MSA files containing all sequences into smaller MSA files and re-writes'
                                             'original MSA file for use with mcorr-pair-sync.')
##Inputs##
parser.add_argument('MSApath', help='path to MSA file to be split')
parser.add_argument('MSAname', help='Name of MSA file to be split')
parser.add_argument('outdir', help='Name of output directory')
parser.add_argument('splits', type = int, help='Number of MSA files to split into')
args = parser.parse_args()
MSA_file = args.MSApath
MSA_name = args.MSAname
outdir = args.outdir
splits = args.splits
###make an output directory
if not os.path.exists(outdir):
    os.makedirs(outdir)

#get the strain set

with open(MSA_file, "r") as MSA:
    strain_set = set()
    for position, ln in enumerate(MSA):
        if ln.startswith(">"):
            strain = str.rstrip(ln.split(' ')[2])
            strain_set.add(strain)

strains = list(strain_set)
print("%s total strains" %str(len(strains)))
##split the set into equal parts
split_strains = np.array_split(strains, splits)
###
print("%s splits performed" %str(len(split_strains)))

##read the master MSA into a file
with open(MSA_file, 'r') as master:
    master_full = master.readlines()

print('re-writing the original MSA ...')
##remake the original MSA file in the same format
with open(MSA_file, "r") as MSA:
    jeans = []
    for position, ln in enumerate(MSA):
        if ln.startswith(">"):
            strain = str.rstrip(ln.split(' ')[2])
            if strain in strains:
                gene = ln.split(' ')[0].split('|')[1]
                jeans.append((position, gene, strain))
outputMSA = os.path.join(outdir, MSA_name)
with open(outputMSA, 'w+') as cluster_MSA:
    lastgene = jeans[0]
    for gene in tqdm(jeans):
        if gene[1] != lastgene[1]:
            cluster_MSA.write('=\n')
            lastgene = gene
        header = master_full[gene[0]]
        seq = master_full[gene[0]+1]
        cluster_MSA.write(header)
        cluster_MSA.write(seq)
    cluster_MSA.write('=')


##now split everything else

print('splitting the MSA ...')

for i in tqdm(np.arange(0, len(split_strains))):
    # cluster = cluster_df[cluster_df['cluster_ID'] == i]
    # clusterstrains = cluster['strain'].tolist()
    clusterstrains = split_strains[i]
    ##get positions of headers, the gene names, and the strain names
    with open(MSA_file, "r") as MSA:
        jeans = []
        for position, ln in enumerate(MSA):
            if ln.startswith(">"):
                strain = str.rstrip(ln.split(' ')[2])
                if strain in clusterstrains:
                    gene = ln.split(' ')[0].split('|')[1]
                    jeans.append((position, gene, strain))
    outputMSA = os.path.join(outdir, MSA_name + '_split' + str(i))
    with open(outputMSA, 'w+') as cluster_MSA:
        lastgene = jeans[0]
        for gene in jeans:
            if gene[1] != lastgene[1]:
                cluster_MSA.write('=\n')
                lastgene = gene
            header = master_full[gene[0]]
            seq = master_full[gene[0]+1]
            cluster_MSA.write(header)
            cluster_MSA.write(seq)
        cluster_MSA.write('=')
