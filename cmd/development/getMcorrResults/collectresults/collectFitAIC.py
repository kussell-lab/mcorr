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
# from mcorr.fit_data import FitDatas
# from mcorr.corr_res import read_corr
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
    parser.add_argument("--corr_dir", default="", help="corr profile directory (if they aren't in the same directory")

    ##define commandline args as variables
    args = parser.parse_args()
    file_dir = args.file_dir
    file_name = args.file_name
    out_dir = args.out_dir
    corr_dir = args.corr_dir
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
    fit_success = []
    #rsq = []

    for csvfile in tqdm(csvpaths):
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
        if "FLEX" in csvfile:
            genome.append("FLEX")
            g = "FLEX"
        else:
            genome.append("CORE")
            g = "CORE"
        ## get the correlation profile to compute rsquared
        # corrcsv = cluster_name + "_" + g + "_XMFA_OUT.csv"
        # if corr_dir == "":
        #     dir = os.getcwd()
        # else:
        #     dir = corr_dir
        # corr_file = os.path.join(dir, cluster_name, corrcsv)
        # # read correlation results
        # corr_results = read_corr(corr_file)
        # fitdatas = FitDatas(corr_results, 3, 300)
        # all = fitdatas.get("all")
        # y = all.yvalues
        # ##get the total variance
        # deltay = y - np.mean(y)
        # SStot = np.sum(deltay**2)
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
        #rsq.append(1-float(params["chisq"])/SStot)
        ##check if the fit succeeded
        stats = open(csvfile)
        for i, line in enumerate(stats):
            if i == 1:
                terms = line.rstrip().split(",")
                fit_success.append(terms[1])

    data = list(zip(cluster_names, cluster_types, genome, modelAIC, lineAIC,
                    d_s, theta_s, fbar, phi_s, theta_p, phi_p, c, d_theta_p,
                    d_theta_s, chisq, red_chisq, fit_success))
    outdat = pd.DataFrame(data, columns=["cluster", "type", "genome", "recombo_AIC", "zero-recombo_AIC",
                                         "ds", "thetaS", "f", "phiS", "thetaP", "phiP", "c", "d_thetaP",
                                         "d_thetaS", "chisq", "red-chisq", "fit_success"])

    now = datetime.datetime.now()
    outpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_'+file_name+'.csv')
    outdat.to_csv(outpath)

if __name__ == "__main__":
    main()