import numpy as np
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
import pandas as pd
from tqdm import tqdm
import os

"""
Goal:

This code is for splitting MSA files containing all sequences in a
WGS dataset into flat sequence clusters using pairwise distances

"""
""""
Inputs
"""
##mcorr-pair output csv
mcp = '/scratch/aps376/APS141_1008_Archive/saureus_mp_MASTER_XMFA_OUT.csv'
##choose what percentile cutoff of the pairwise distances
# you'll use to make flat clusters
cutoff = 10
#Where the MSA master file is
##MAKE SURE TO USE THE SORTED ONE!!
MSA_file = '/scratch/aps376/APS141_1008_Archive/APS141_saureus_MASTER_MSA'

###output directory if you so please (but you should please)

outdir = '/scratch/aps376/APS144_1008_saureus_Archive/cutoff_10pt/'
##save directory for the distance matrices
archive = '/scratch/aps376/APS144_1008_saureus_Archive/'

"""
function for converting mcorr-pair output to a distance matrix
returns a list of strain names and the distance matrix
"""
def mptodm(csv):
    #returns two distance matrices from an mcorr-pair .csv output file
    ##one is a distance matrix using pairwise diversity or d_sample
    #called dm_d. The other is a distance matrix using theta_pool
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    ##read in mcorr-pair output
    dat = pd.read_csv(csv)
    dat = dat[dat['b']!= 'all']
    #first make full matrices of both parameters that are
    #symmetric about the diagonal (which is 0)
    #get the names of each strain in the rows (which are identical to collumn labels)
    pairs = dat['b']
    firstpair = pairs[0]
    pairname1 = firstpair.split("_vs_")
    firstrowname = pairname1[0]
    firstrow = dat[dat['b'].str.contains(firstrowname)]
    firstrow = firstrow.reset_index(drop = True)
    rownames = []
    rownames.append(firstrowname)
    for i in np.arange(0, len(firstrow)):
        pair = firstrow['b'][i]
        pairname = pair.split("_vs_")
        if pairname[0] != firstrowname:
            rownames.append(pairname[0])
        else:
            rownames.append(pairname[1])
    ##get d_sample
    firstrow_d =  np.append([0.0], firstrow['m'], axis = 0)
    #store the names
    names = rownames
    ##store the first row of the distance matrix
    fullm_d = firstrow_d

    ##now get the rest of the rows
    #j = 1
    for j in tqdm(np.arange(1,len(rownames))):
        row = dat[dat['b'].str.contains(rownames[j])]
        row = row.reset_index(drop = True)
        ##second row first element
        row_element1 = row[row['b'].str.contains(rownames[0])]
        row_d = row_element1['m']
        for rowname in rownames:
            if rowname == rownames[0]:
                continue
            if rowname == rownames[j]:
                row_d = np.append(row_d, 0.0)
            else:
                row_element = row[row['b'].str.contains(rowname)]
                row_d = np.append(row_d, row_element['m'])
        fullm_d = np.row_stack((fullm_d, row_d))


    ##return dm_d and dm_theta
    return names, fullm_d

### if the distance matrix doesn't exist, make it
if not os.path.exists(archive+'fullm_d'):
    ##get the distance matrix
    names, fullm_d = mptodm(mcp)
    ### save the distance matrices for later (could be useful for future analysis)
    np.save(archive+'fullm_d', fullm_d)
    np.save(archive+'names', names)

###load up the distance matrix and the names
fullm_d = np.load(archive+'fullm_d.npy')
names = np.load(archive+'names.npy')
"""
convert distance matrices to flat dist matrices
(returns flat dist matrix)
"""
def dm2flatdm(fullm_d):
    import numpy as np
    from tqdm import tqdm
    ##make flat distance matrices for scipy
    flatdm_d = []
    for i in tqdm(np.arange(0, len(fullm_d))):
        flatrow_d = fullm_d[i][i+1:].tolist()
        flatdm_d = flatdm_d+flatrow_d
    return flatdm_d

flatdm_d = dm2flatdm(fullm_d)

###make your clusters
## get some stats, generate a submission list
flat = np.asarray(flatdm_d)
pairwise = flat[flat != 0]

Zavg_d = linkage(flatdm_d, 'average')
c, coph_dists = cophenet(Zavg_d, flatdm_d)
print('Cophenetic: '+str(c))

clusters = fcluster(Zavg_d, np.percentile(pairwise, cutoff), criterion='distance')
unique, counts = np.unique(clusters, return_counts=True)
clusterlist = dict(zip(unique, counts))

totclusters = 0
allstrains = 0
analyzedstrains = 0
submissionlist = []
for i in unique:
    allstrains = allstrains+clusterlist[i]
    if clusterlist[i] == 1:
        continue
    totclusters = totclusters + 1
    analyzedstrains = analyzedstrains+clusterlist[i]
    submissionlist.append(i)
###print how many clusters and strains will be analyzed here
###and a submission list for the clusters (you'll use this later with mcorr)
print('Clusters with more than one sequence: '+str(totclusters))
print('All strains: '+str(allstrains))
print('Analyzed strains: '+str(analyzedstrains))
print('Submission list:')
print(submissionlist)
print('Submission list sans commas:')
print(*submissionlist, sep=' ')

"""
Make MSA files for each cluster
"""

strains = names

cluster_df = pd.DataFrame(list(zip(clusters.tolist(), strains.tolist())),
                          columns = ['cluster_ID', 'strain'])

with open(MSA_file, 'r') as master:
    master_full = master.readlines()

if not os.path.exists(outdir):
    os.makedirs(outdir)
for i in tqdm(set(clusters)):
    cluster = cluster_df[cluster_df['cluster_ID'] == i]
    clusterstrains = cluster['strain'].tolist()
    if len(clusterstrains) < 2:
        continue
    ##get positions of headers, the gene names, and the strain names
    with open(MSA_file, "r") as MSA:
        jeans = []
        for position, ln in enumerate(MSA):
            if ln.startswith(">"):
                strain = str.rstrip(ln.split(' ')[2])
                if strain in clusterstrains:
                    gene = ln.split(' ')[0].split('|')[1]
                    jeans.append((position, gene, strain))

    with open(outdir+'MSA_cluster'+str(i), 'w+') as cluster_MSA:
        lastgene = jeans[0]
        for gene in jeans:
            if gene[1] != lastgene[1]:
                cluster_MSA.write('=\n')
                lastgene = gene
            header = master_full[gene[0]]
            seq = master_full[gene[0]+1]
            cluster_MSA.write(header)
            cluster_MSA.write(seq)

""""
A function for splitting cluster MSAs into core and flex MSA files
"""
def coreflexsplit(inputdir, outdir, threshold):
    ##(1) input folder with the MSA files
    ##(2) output folder for split files
    ##(3) cutoff to be considered a core or flexible gene

    import os
    import glob
    import sys
    import pandas as pd
    from tqdm import tqdm
    ##change directories to the input directory
    os.chdir(inputdir)
    ##if the output folder doesn't exist, make it
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    #Create a dataframe called 'get_mapped' that contains the percentage of
    #sequence cluster samples that map to each gene
    get_mapped = pd.DataFrame()
    ###### read through all of the MSA files in the folder
    for refout in tqdm(glob.glob("MSA*")):
        msaname = refout.split('_')
        #Make a list for core and a list for flex genes
        with open(refout,'r') as cluster_file:
            cluster_lines = cluster_file.readlines()

        #Gene information
        gene_end = [i for i, x in enumerate(cluster_lines) if x == '=\n']
        gene_start = [-1]+gene_end[:-1]
        gene_name = [cluster_lines[name_ind-2].replace('|',' ').split()[1] for name_ind in gene_end]

        for g, name in enumerate(gene_name):
            gene_range = cluster_lines[gene_start[g]+1:gene_end[g]+1]
            seq_lines = [len(set(x)) > 2 for i, x in enumerate(gene_range) if x[0] != '>']
            ## looks like should be len(seq_lines)-1 in the bottom
            mapped_percentage = sum(seq_lines)*1/(len(seq_lines)-1)
            get_mapped.at[name, msaname[1]] = mapped_percentage

    #Inelegant. Verging on trash... but it works!
    #This snippet determines what genes are CORE genes by seeing what genes have
    #enough mapped samples to be above the threshold. It then creates a DataFrame with
    #a 'True' value when that gene is represented across all sequence clusters.
    core_genes = (get_mapped>threshold)*1
    core_genes = core_genes.sum(axis=1).astype(int)
    core_genes = core_genes == len(get_mapped.columns)

    #Export get_mapped
    #This exported dataframe has a gene name index and a column for each sequence cluster.
    #The values are what percentage of samples, for that sequence cluster, had that gene.
    get_mapped.to_csv(outdir+'gene_percentage_all.csv')

    #Create CORE and FLEX files for each sequence cluster
    for refout in tqdm(glob.glob("MSA*")):

        #Make a list for core and a list for flex genes
        with open(refout,'r') as cluster_file:
            cluster_lines = cluster_file.readlines()

        #Gene core vs flex
        gene_end = [i for i, x in enumerate(cluster_lines) if x == '=\n']
        gene_start = [-1]+gene_end[:-1]
        gene_name = [cluster_lines[name_ind-2].replace('|',' ').split()[1] for name_ind in gene_end]

        core_lines = []
        flex_lines = []
        for g, name in enumerate(gene_name):
            gene_range = cluster_lines[gene_start[g]+1:gene_end[g]+1]
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
            ##usually ref_dir, changed for this
            msaname = refout.split('_')
            refout_CF = outdir+msaname[0]+'_'+CF_names[CF]+'_'+msaname[1]
            with open(refout_CF, 'w+') as REFGEN_out:
                for line in CF_lists[CF][:]:
                    REFGEN_out.write(line)

inputdir = outdir
threshold = 0.9

coreflexsplit(inputdir, outdir, threshold)


