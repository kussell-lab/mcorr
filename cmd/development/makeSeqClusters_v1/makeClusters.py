#!/usr/bin/env python3
import numpy as np
import argparse
from scipy.cluster.hierarchy import fcluster, linkage, cophenet
import os
import pandas as pd
from scipy.spatial.distance import squareform

def main():

    parser = argparse.ArgumentParser(description="Takes distance matrix and strain list and breaks them into clusters of strains"
                                                 +"using average linkage")

    ##distance matrix
    parser.add_argument("distance_matrix", help="distance matrix from mcorr-pair as a npy matrix array")
    ##strain names
    parser.add_argument("strains", help="strain names ordered to match distance matrix")
    ##percentile cutoff for pairwise distances
    parser.add_argument("cutoff", type=int, help="cutoff percentile of distances used for making clusters (%)")
    ##define commandline args as variables
    args = parser.parse_args()
    fullm_d = np.load(args.distance_matrix,
                      allow_pickle=True,
                      fix_imports=True)
    cutoff = args.cutoff
    names = np.loadtxt(args.strains, dtype=str)
    archive = os.getcwd()
    """""
    inputs, soon to be converted to argparse format
    """
    # archive = '/Volumes/GoogleDrive/My Drive/hpc/recombo/APS150_SP_distances/'
    # fullm_d = np.load(args.distance_matrix,
    #                   allow_pickle=True,
    #                   fix_imports=True)
    # names = np.loadtxt(archive+'strains', dtype=str)
    # MSA_file = '/Volumes/GoogleDrive/My Drive/hpc/recombo/APS150_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER'
    # cutoff = 10

    """
    convert distance matrices to flat dist matrices
    (returns flat dist matrix)
    """
    # ##get indices of upper triangle matrix
    # iu = np.triu_indices(len(fullm_d), k=1)
    # ##get a "condensed distance matrix" for scipy
    # flatdm_d = fullm_d[iu]

    ##convert to condensed distance matrix
    flatdm_d = squareform(fullm_d)
    #print(flatdm_d[-5:])

    ##make the clusters, get some stats, and generate a submission list for the cluster
    ##cluster via average linkage
    Zavg_d = linkage(flatdm_d, 'average')
    ##get the cophenetic
    c, coph_dists = cophenet(Zavg_d, flatdm_d)

    clusters = fcluster(Zavg_d, np.percentile(flatdm_d, cutoff), criterion='distance')
    unique, counts = np.unique(clusters, return_counts=True)
    clusterlist = dict(zip(unique, counts))

    ###print how many clusters and strains will be analyzed here
    ###and a submission list for the clusters (you'll use this later with mcorr)
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
    stats = os.path.join(archive, "stats.txt")
    f = open(stats, "w+")
    f.write('Cophenetic: '+str(c)+'\n')
    f.write('Clusters with more than one sequence: '+str(totclusters)+'\n')
    f.write('All strains: '+str(allstrains)+'\n')
    f.write('Analyzed strains: '+str(analyzedstrains)+'\n')
    f.write('Submission list:\n')
    s = ", ".join(map(str, submissionlist))
    f.write(s+'\n')
    f.write('Submission list sans commas:\n')
    s = " ".join(map(str, submissionlist))
    f.write(s)
    f.close()

    cluster_df = pd.DataFrame(list(zip(clusters, names)),
                              columns = ['cluster_ID', 'strain'])
    grouped = cluster_df.groupby("cluster_ID")

    filtered_df = grouped.filter(lambda x: len(x) > 1.)
    filtered_df = filtered_df.reset_index(drop=True)
    filtered_df = filtered_df[['strain', 'cluster_ID']]
    # clusterdict = dict(zip(names, clusters))
    # clusterlist = pd.DataFrame(list(zip(names, clusters)))
    # saveclusters = os.path.join(archive, "cluster_dict.npy")
    # np.save(saveclusters, clusterdict)
    saveclusters = os.path.join(archive, "cluster_list")
    filtered_df.to_csv(path_or_buf=saveclusters, sep=",", header=False, index=False)

if __name__ == "__main__":
    main()

# ##get positions of headers, the gene names, and the strain names
# with open(MSA_file, "r") as MSA:
#     jeans = []
#     for position, ln in enumerate(MSA):
#         if ln.startswith(">"):
#             strain = str.rstrip(ln.split(' ')[2])
#             gene = ln.split(' ')[0].split('|')[1]
#             jeans.append((position, gene, strain, clusterdict[strain]))
#
# jeans_df = pd.Dataframe(jeans, columns=['position', 'gene', 'strain', 'cluster'])
# print('mammajammalamma')

##store the lines of the master MSA as a list
# with open(MSA_file, 'r') as master:
#     master_full = master.readlines()
# print("going to start daskin")
# start_time = time.time()
# master_full = pd.read_csv(MSA_file, header=None, sep="\n",dtype=str)
# print("--- %s seconds ---" % (time.time() - start_time))
# headers = master_full[master_full.str.contains(">")]
# jeans = headers.lines.str.split(expand=True)
# jeans.rename(columns={"0": "gene", "1": "pos", "2": "strain"})
# seqs = np.array(jeans["strain"])
# straincluster = np.empty(len(jeans), dtype= int)
# for i in np.arange(0, len(seqs)):
#     straincluster[i] = clusterdict[seqs[i]]
# jeans["cluster"] = straincluster
# print("jamma committee")