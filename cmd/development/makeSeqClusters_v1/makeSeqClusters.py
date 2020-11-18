#!/usr/bin/env python3
import time
import timeit
import numpy as np
import argparse
import os
"""
this commandline program controls everything else
"""

parser = argparse.ArgumentParser(description="Takes the output of mcorr-pair, clusters sequences, and breaks clusters into core and flexible genomes")
parser.add_argument("mcorr_pair_csv", help="the output of pairwise distances from mcorr-pair")
parser.add_argument("working_dir", help="the working space and output directory")
parser.add_argument("cutoff", type=int, help="cutoff percentile of pairwise distances to use for making flat clusters (%)")
parser.add_argument("master_msa", type=str, help="master msa file for all strain sequences")
parser.add_argument("threshold", type=int, help="threshold percentage above which you're considered a core gene (%)")
##define commandline args as variables
args = parser.parse_args()
mcp = args.mcorr_pair_csv
wrkdir = args.working_dir
percentile = args.cutoff
msa = args.master_msa
t = args.threshold


# fullm_d = np.load(args.distance_matrix,
#                   allow_pickle=True,
#                   fix_imports=True)
# cutoff = args.cutoff
# names = np.loadtxt(args.strains, dtype=str)
# archive = os.getcwd()
start_time = time.time()

os.chdir(wrkdir)

#Step 1: take mcorr-pair output and translate to distance matrix and a list of strains
print("converting mcorr-pair to matrix ...")
os.system("mcorr-dm %s distancematrix" %mcp)

#Step 2: take the distance matrix and list of names
print("clustering sequences")
os.system("python3 makeSeqClusters.py distancematrix strains %s" %percentile)

#Step 3: make sequence cluster MSA files and core and flexible MSA files
print("writing cluster MSA files ...")
os.system("writeclusters $s cluster_list %s" % (msa, t))
print("Done with making clusters, time to boogie")
print("--- %s minutes ---" % (time.time() - start_time)/60)