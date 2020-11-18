#!/usr/bin/env python3
import time
import argparse
import os
from . import clusterSequences as cS
"""
this commandline program controls everything else
"""
def main():
    parser = argparse.ArgumentParser(description="Takes the output of mcorr-pair, clusters sequences,\
                                                    and breaks clusters into core and flexible genomes")
    parser.add_argument("mcorr_pair_csv", help="the output of pairwise distances from mcorr-pair")
    parser.add_argument("working_dir", help="the working space and output directory")
    parser.add_argument("--percentile", type=float, default=10, help="cutoff percentile of pairwise distances to use for making flat clusters (Default: 10)")
    parser.add_argument("master_msa", type=str, help="master msa file for all strain sequences")
    parser.add_argument("--threshold", type=int, default=90, help="threshold percentage above which you're considered a core gene (Default: 90)")
    ##define commandline args as variables
    args = parser.parse_args()
    mcp = args.mcorr_pair_csv
    wrkdir = args.working_dir
    percentile = args.percentile
    msa = args.master_msa
    t = args.threshold

    start_time = time.time()

    os.chdir(wrkdir)

    #Step 1: take mcorr-pair output and translate to distance matrix and a list of strains
    print("converting mcorr-pair to matrix ...")
    os.system("mcorr-dm %s distancematrix" %mcp)

    #Step 2: take the distance matrix and list of names
    print("clustering sequences ...")
    cS.clusterSequences("distancematrix.npy", "strains", percentile)

    #Step 3: make sequence cluster MSA files and core and flexible MSA files
    print("writing cluster MSA files ...")
    os.system("write-cluster-msa %s cluster_list --core-cutoff=%s" % (msa, t))
    print("Done with making clusters, time to boogie")
    print("Total run time: %s minutes" % str((time.time() - start_time)/60))

if __name__ == "__main__":
    main()