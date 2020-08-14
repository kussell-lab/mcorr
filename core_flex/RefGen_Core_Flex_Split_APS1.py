import os
import sys
import pandas as pd

###
#Background
#This script assumes you've already run the 'pipeline' bash script
#After you've run the 'pipeline' script you should have a series of folders for each serotype
#This will crawl through those folders and create CORE and FLEX gene XMFA files within them

#Notes
#This program has room for improvement. It works, but it's a little disorganized.
###

#INPUTS
#sero_filelist is a text file with each line giving the name of a serotype
sero_filelist = '/scratch/aps376/sero_list_8'
#ref_dir is the directory that contains the all of the serotype folders generated by pipeline.sh
ref_dir = '/scratch/aps376/Archive/'
#The threshold defining what percentage of samples need to map before a gene is considered part of the serotype
threshold = 0.9

#Get the list of serotypes
with open(sero_filelist,'r') as f:
    seros = f.read().splitlines()

#Create a dataframe called 'get_mapped' that contains the percentage of
#serotype samples that map to each gene
get_mapped = pd.DataFrame()
for sero_name in seros:
    #Where the REFGEN sero file is
    refout = ref_dir + sero_name + '_OUT/REFGEN_' + sero_name

    #Make a list for core and a list for flex genes
    with open(refout,'r') as sero_file:
        sero_lines = sero_file.readlines()

    #Gene information
    gene_end = [i for i, x in enumerate(sero_lines) if x == '=\n']
    gene_start = [-1]+gene_end[:-1]
    gene_name = [sero_lines[name_ind-2].replace('|',' ').split()[1] for name_ind in gene_end]
    #sero_df = pd.DataFrame(columns=[sero_name],index=gene_name)

    for g, name in enumerate(gene_name):
        gene_range = sero_lines[gene_start[g]+1:gene_end[g]+1]
        seq_lines = [len(set(x)) > 2 for i, x in enumerate(gene_range) if x[0] != '>']
        mapped_percentage = sum(seq_lines)*1/len(seq_lines)
        get_mapped.at[name,sero_name] = mapped_percentage

#Inelegant. Verging on trash... but it works!
#This snippet determines what genes are CORE genes by seeing what genes have
#enough mapped samples to be above the threshold. It then creates a DataFrame with
#a 'True' value when that gene is represented across all serotypes.
core_genes = (get_mapped>threshold)*1
core_genes = core_genes.sum(axis=1).astype(int)
core_genes = core_genes == len(get_mapped.columns)

#Export get_mapped
#This exported dataframe has a gene name index and a column for each serotype.
#The values are what percentage of samples, for that serotype, had that gene.
get_mapped.to_csv('gene_percentage_all.csv', sep='\t')

#Create CORE and FLEX files for each serotype
for sero_name in seros:
    #Where the REFGEN sero file is
    refout = ref_dir + sero_name + '_OUT/REFGEN_' + sero_name

    #Make a list for core and a list for flex genes
    with open(refout,'r') as sero_file:
        sero_lines = sero_file.readlines()

    #Gene core vs flex
    gene_end = [i for i, x in enumerate(sero_lines) if x == '=\n']
    gene_start = [-1]+gene_end[:-1]
    gene_name = [sero_lines[name_ind-2].replace('|',' ').split()[1] for name_ind in gene_end]

    core_lines = []
    flex_lines = []
    for g, name in enumerate(gene_name):
        gene_range = sero_lines[gene_start[g]+1:gene_end[g]+1]
        seq_lines = [len(set(x)) > 2 for i, x in enumerate(gene_range) if x[0] != '>']
        #mapped_percentage = sum(seq_lines)*1/len(seq_lines)
        if core_genes[name]:
            #CORE gene
            core_lines.extend(gene_range)
        else:
            #FLEX gene
            flex_lines.extend(gene_range)

    CF_lists = [core_lines,flex_lines]
    for CF in range(0,2):
        CF_names = ['CORE','FLEX']
        #Where to write files to
        refout_CF = ref_dir + sero_name + '_OUT/REFGEN_' + CF_names[CF] + '_' + sero_name
        with open(refout_CF, 'a+') as REFGEN_out:
            for line in CF_lists[CF][:]:
                REFGEN_out.write(line)

print('done making banana splits')