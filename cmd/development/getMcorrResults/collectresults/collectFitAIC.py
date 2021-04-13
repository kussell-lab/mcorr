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

    ##make the out directory if it doesn't exist
    mkdir_p(out_dir)

    ##get cluster names from the stats file
    ##remove this if you want to do a local test
    globpath = os.path.join(file_dir, "cluster*", "*"+suffix)
    csvpaths = glob.glob(globpath)
    cluster_names = []
    cluster_types = []
    genome = []
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

    for csvfile in tqdm(csvpaths):
        terms = csvfile.split("/")
        cluster_name = terms[len(terms)-2]
        cluster_names.append(cluster_name)
        terms = cluster_name.split("_")
        if len(terms) > 1:
            cluster_types.append("between")
        else:
            cluster_types.append("within")

        dat = pd.read_csv(csvfile, skiprows=10)
        AIC = np.array(dat["AIC"])
        modelAIC.append(AIC[0])
        lineAIC.append(AIC[1])
        if "FLEX" in csvfile:
            genome.append("FLEX")
        else:
            genome.append("CORE")
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

    data = list(zip(cluster_names, cluster_types, genome, modelAIC, lineAIC,
                    d_s, theta_s, fbar, phi_s, theta_p, phi_p, c, d_theta_p,
                    d_theta_s, chisq, red_chisq))
    outdat = pd.DataFrame(data, columns=["cluster", "type", "genome", "recombo_AIC", "zero-recombo_AIC",
                                         "d_s", "theta_s", "f", "phi_s", "theta_p", "phi_p", "c", "d_theta_p",
                                         "d_theta_s", "chisq", "red-chisq"])

    now = datetime.datetime.now()
    outpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_'+file_name+'.csv')
    outdat.to_csv(outpath)

if __name__ == "__main__":
    main()