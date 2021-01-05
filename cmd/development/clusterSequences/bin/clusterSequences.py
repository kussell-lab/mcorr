#!/usr/bin/env python3
import time
import numpy as np
import scipy
import argparse
from scipy.cluster.hierarchy import fcluster, linkage, cophenet, dendrogram
import os
import pandas as pd
from scipy.spatial.distance import squareform
import matplotlib as mpl
mpl.use('Agg')
from matplotlib import pyplot as plt
import seaborn as sns

def main():
    parser = argparse.ArgumentParser(description="Takes a distance matrix, clusters sequences, and outputs a cluster list\
                                                    for use with write-cluster-msa ")
    parser.add_argument("distance_matrix", help="symmetric distance matrix as .npy file")
    parser.add_argument("strains", help="list of strains")
    parser.add_argument("dendrogram", help="prefix for dendogram output file")
    parser.add_argument("--percentile", type=float, default=10, help="cutoff percentile of pairwise distances to use for making flat clusters (Default: 10)")
    ##define commandline args as variables
    args = parser.parse_args()
    distance_matrix = args.distance_matrix
    strains = args.strains
    cutoff = args.percentile
    outfile = args.dendrogram
    ####
    start_time = time.time()
    fullm_d = np.load(distance_matrix,
                      allow_pickle=True,
                      fix_imports=True)
    #cutoff = args.cutoff
    names = np.loadtxt(strains, dtype=str)
    archive = os.getcwd()
    cutoff = float(cutoff)

    ##make sure your matrix is symmetric
    if (fullm_d.transpose() == fullm_d).all():
        print("It's symmetric")
    else:
        print("Not symmetric")

    ##convert to condensed distance matrix
    flatdm_d = squareform(fullm_d)

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
    stats = os.path.join(archive, str(int(cutoff))+"th_percentile_stats.txt")
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
    saveclusters = os.path.join(archive, "cluster_list")
    filtered_df.to_csv(path_or_buf=saveclusters, sep=",", header=False, index=False)
    ##get colors for the dendrogram
    colors = sns.color_palette("nipy_spectral", totclusters)
    accents = []
    for color in colors:
        this = mpl.colors.to_hex(color)
        accents.append(this)
    scipy.cluster.hierarchy.set_link_color_palette(accents)
    ##make a dendrogram plot as a pdf
    percentage = "{:.4f}".format(np.percentile(flatdm_d, cutoff))
    mpl.rcParams['lines.linewidth'] = 0.5
    fig, ax = plt.subplots(figsize=(10, 30))
    title = "%sth percentile cutoff ($P_{10}=%s$)" % (str(int(cutoff)), percentage)
    plt.title(title, fontsize=16)
    plt.ylabel('')
    plt.xlabel('distance')
    dendrogram(
        Zavg_d,
        no_labels=True,
        color_threshold=np.percentile(flatdm_d, cutoff),
        orientation='left',
    )

    plt.axvline(x=np.percentile(flatdm_d, cutoff), c='k', alpha=0.7, linestyle='--')

    plt.savefig(outfile+'.pdf')

    print("Time to cluster sequences: %s s" % str(time.time() - start_time))

if __name__ == "__main__":
    main()