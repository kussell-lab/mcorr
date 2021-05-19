#!/usr/bin/env python3
import datetime
import glob
import time
import argparse
import os
import csv
from tqdm import tqdm
import pandas as pd
from itertools import combinations
import numpy as np
import linecache
#from .functions import mkdir_p
"""
Commandline program for collecting the model comparison results from mcorrFitCompare
"""

def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)

def main():
    parser = argparse.ArgumentParser(description="Collect results from fitting and the model comparison from mcorrFitCompare")
    parser.add_argument("--file_dir", default="current", help="specify directory with output from mcorr if not current directory")
    parser.add_argument("--out_dir", default="current", help="optional out directory for output csv")
    parser.add_argument("--suffix", default="", help="if the goodness-of-fit stats sheet has a different suffix")
    parser.add_argument("file_name", help="prefix for resultant csv file (timestamp will be included)")

    ##define commandline args as variables
    args = parser.parse_args()
    file_dir = args.file_dir
    file_name = args.file_name
    out_dir = args.out_dir
    suffix = args.suffix

    if suffix == "":
        suffix = '_comparemodels.csv'

    ##for local testing
    # file_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # out_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # file_name = 'mcorr_res_test'
    # suffix = '_FIT_OUT_lmfit_report.txt'
    # clusters = [8, 221]

    if file_dir == "current":
        file_dir = os.getcwd()
    if out_dir == "current":
        out_dir = os.getcwd()

    ##make a bin list
    globpath = os.path.join(file_dir, "bin*")
    binpaths = glob.glob(globpath)
    all_bins = []
    for bin in binpaths:
        terms = bin.split("/")
        bindir = terms[len(terms)-1]
        bin = bindir.split("_")[1]
        all_bins.append(bin)
    ##make the out directory if it doesn't exist
    mkdir_p(out_dir)

    ##check which outputs are missing and print a list
    print("make cluster list ...")
    globpath = os.path.join(file_dir, "bin_"+all_bins[0], "cluster*")
    clusterpaths = glob.glob(globpath)
    all_clusters = []
    ##also make a set of bins
    for cluster in clusterpaths:
        terms = cluster.split("/")
        cluster_name = terms[len(terms)-1]
        all_clusters.append(cluster_name)
    allclusters = set(all_clusters)

    ##get cluster names from the stats file and to check which outputs are missing
    ##remove this if you want to do a local test
    missing_clusters = []
    missing_bins = []
    for bin in all_bins:
        globpath = os.path.join(file_dir, "bin_" + bin, "cluster*", "*"+suffix)
        csvpaths = glob.glob(globpath)
        outputs = []
        for cluster in csvpaths:
            terms = cluster.split("/")
            cluster_name = terms[len(terms)-2]
            outputs.append(cluster_name)
        outs = set(outputs)
        #### see what's missing ...
        missing = allclusters.difference(outs)
        for cluster in missing:
            missing_clusters.append(cluster)
            missing_bins.append(bin)
    ##print it to a list so we can re-run these jobs
    data = list(zip(missing_clusters, missing_bins))
    outdat = pd.DataFrame(data, columns=["cluster", "bin"])
    ##get the time
    now = datetime.datetime.now()
    outpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_unfinished_jobs.csv')
    outdat.to_csv(outpath)
    ########
    cluster_names = []
    cluster_types = []
    modelAIC = []
    lineAIC = []
    d_s = []
    theta_s = []
    fbar = []
    phi_s = []
    theta_p = []
    phi_p = []
    c = []
    d_theta_p = []
    d_theta_s = []
    chisq = []
    red_chisq = []
    fit_success = []
    binrange = []

    for bin in all_bins:
        globpath = os.path.join(file_dir, "bin_" + bin, "cluster*", "*"+suffix)
        csvpaths = glob.glob(globpath)
        for csvfile in tqdm(csvpaths):
            binrange.append(bin)
            terms = csvfile.split("/")
            cluster_name = terms[len(terms)-2]
            cluster_names.append(cluster_name)
            terms = cluster_name.split("_")
            if len(terms) > 1:
                cluster_types.append("Btwn clusters")
            else:
                cluster_types.append("W/n cluster")

            dat = pd.read_csv(csvfile, skiprows=10)
            AIC = np.array(dat["AIC"])
            modelAIC.append(AIC[0])
            lineAIC.append(AIC[1])
            ##collect params
            params = dat[dat["recombination"]=="recombo"]
            d_s.append(float(params["d_s"]))
            theta_s.append(float(params["theta_s"]))
            fbar.append(float(params["f"]))
            phi_s.append(float(params["phi_s"]))
            theta_p.append(float(params["theta_p"]))
            phi_p.append(float(params["phi_p"]))
            c.append(float(params["c"]))
            d_theta_p.append(float(params["d_theta_p"]))
            d_theta_s.append(float(params["d_theta_s"]))
            chisq.append(float(params["chisq"]))
            red_chisq.append(float(params["red-chisq"]))
            ##check if the fit succeeded
            stats = open(csvfile)
            for i, line in enumerate(stats):
                if i == 1:
                    terms = line.rstrip().split(",")
                    fit_success.append(terms[1])

    data = list(zip(cluster_names, cluster_types, binrange, modelAIC, lineAIC,
                    d_s, theta_s, fbar, phi_s, theta_p, phi_p, c, d_theta_p,
                    d_theta_s, chisq, red_chisq, fit_success))
    outdat = pd.DataFrame(data, columns=["cluster", "type", "bin", "recombo_AIC", "zero-recombo_AIC",
                                         "ds", "thetaS", "f", "phiS", "thetaP", "phiP", "c", "d_thetaP",
                                         "d_thetaS", "chisq", "red-chisq", "fit_success"])

    outpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_'+file_name+'.csv')
    outdat.to_csv(outpath)

if __name__ == "__main__":
    main()