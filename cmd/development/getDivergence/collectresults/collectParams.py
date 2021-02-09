#!/usr/bin/env python3
import datetime
import time
import argparse
import os
import csv
from tqdm import tqdm
import pandas as pd
from itertools import combinations
import numpy as np
import linecache
from collectresults.functions import mkdir_p, get_fitstats, get_numalns
"""
Commandline program for collecting divergence and bootstrap interval outputs from mcorr-fit
"""

def main():
    parser = argparse.ArgumentParser(description="Collect all inferred params for sequence clusters from mcorr-fit\n"+
                                     "outputs are (1) csv file params, (2) a list of clusters with incomplete\n"+
                                                 " jobs named 'DATETIMESTAMP_incomplete.txt'")
    parser.add_argument("--file_dir", default="current", help="specify directory with output from mcorr if not current directory")
    parser.add_argument("--out_dir", default="current", help="optional output directory for divergence csv")
    parser.add_argument("--suffix", default="", help="if the lmfit goodness-of-fit stats sheet has a different suffix")
    parser.add_argument("stats_sheet", help=" name of stats csv file from clusterSequences.py which has the list of clusters")
    parser.add_argument("file_name", help="prefix for resultant csv file (timestamp will be included)")

    ##define commandline args as variables
    args = parser.parse_args()
    file_dir = args.file_dir
    file_name = args.file_name
    stats_sheet = args.stats_sheet
    out_dir = args.out_dir
    suffix = args.suffix

    ##for local testing
    # file_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # out_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # file_name = '0122_mcorr_res_test'
    # suffix = '_FIT_OUT_lmfit_report.txt'
    # clusters = [8, 221]
    #stats_sheet = "20th_percentile_>=10strains_stats.csv"

    if file_dir == "current":
        file_dir = os.getcwd()
    if out_dir == "current":
        out_dir = os.getcwd()

    ##make the out directory if it doesn't exist
    mkdir_p(out_dir)

    ##get cluster names from the stats file
    ##remove this if you want to do a local test
    stats = os.path.join(file_dir, stats_sheet)
    with open(stats) as csvDataFile:
        statsdat=[row for row in csv.reader(csvDataFile)]
    clusterlist = statsdat[5][1]
    clusters = [int(s) for s in clusterlist.split(',')]



    ##i am going to refer to these as serotypes bc i'm too lazy to edit this, but know ...
    ##they are sequence clusters!
    serocombs = combinations(clusters, 2)
    combolist = []
    for c in serocombs:
        combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))
    sero_list = []
    for cluster in clusters:
        sero_list.append('cluster'+str(cluster))

    i = 0
    ##check if any runs didn't finish
    repeats = []
    ##ones that weren't fit by curve-fitting
    misfits = []
    ##for just within sero data
    data = []
    for sero in tqdm(sero_list):
        core = True
        flex = True
        core_file = os.path.join(file_dir, str(sero), str(sero)+'_CORE_FIT_OUT_fit_results.csv')
        flex_file = os.path.join(file_dir, str(sero), str(sero)+'_FLEX_FIT_OUT_fit_results.csv')
        corestats = os.path.join(file_dir, str(sero), str(sero)+'_CORE_FIT_OUT_lmfit_report.txt')
        flexstats = os.path.join(file_dir, str(sero), str(sero)+'_FLEX_FIT_OUT_lmfit_report.txt')
        if not os.path.exists(core_file):
            repeats.append((c[0], c[1], 'CORE'))
            core = False
        if not os.path.exists(corestats):
            misfits.append((c[0], c[1], 'CORE'))
            core = False
        if not os.path.exists(flex_file):
            repeats.append((c[0], c[1], 'FLEX'))
            flex = False
        if not os.path.exists(flexstats):
            misfits.append((c[0], c[1], 'FLEX'))
            flex = False

        ##core
        if core:
            core_row = []
            core = pd.read_csv(core_file)
            core_row.append(str(sero))
            core_row.append("CORE")
            core_row.append("within cluster")
            core_data = core[core["group"] == "all"].values.flatten().tolist()
            core_row = core_row + core_data

            ##get fitstats
            points, vars, chi, redchi = get_fitstats(file_dir, sero, "", "CORE", suffix)
            core_row = core_row + [points, vars, chi, redchi]
            data.append(core_row)
            i = i + 1

        ##flex
        if flex:
            flex_row = []
            flex = pd.read_csv(flex_file)
            flex_row.append(str(sero))
            flex_row.append("flex")
            flex_row.append("within cluster")
            flex_data = flex[flex["group"] == "all"].values.flatten().tolist()
            flex_row = flex_row + flex_data
            ##get fitstats
            points, vars, chi, redchi = get_fitstats(file_dir, sero, "", "FLEX", suffix)
            flex_row = flex_row + [points, vars, chi, redchi]
            data.append(flex_row)
            i = i + 1

    if i != 0:
        within = pd.DataFrame(data,
                              columns=['cluster', 'genome', 'type', 'group','d_sample','theta_pool',
                                       'phi_pool','ratio','fbar','c','d_pool',
                                       'd_clonal','theta_s','phi_s', 'bp_analyzed', 'variables', 'chisq', 'red-chisq'])


    "for between sero data"
    data = []

    j = 0

    for c in tqdm(combolist):
        core = True
        flex = True
        core_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE_FIT_OUT_fit_results.csv')
        flex_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX_FIT_OUT_fit_results.csv')
        corestats = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE_FIT_OUT_lmfit_report.txt')
        flexstats = os.path.join(file_dir,  c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX_FIT_OUT_lmfit_report.txt')
        if not os.path.exists(core_file):
            repeats.append((c[0], c[1], 'CORE'))
            core = False
        if not os.path.exists(corestats):
            misfits.append((c[0], c[1], 'CORE'))
            core = False
        if not os.path.exists(flex_file):
            repeats.append((c[0], c[1], 'FLEX'))
            flex = False
        if not os.path.exists(flexstats):
            misfits.append((c[0], c[1], 'FLEX'))
            flex = False
        ##core
        if core:
            core_row = []
            core = pd.read_csv(core_file)
            core_row.append(c[0]+'/'+c[1])
            core_row.append("CORE")
            core_row.append("between cluster")
            core_data = core[core["group"] == "all"].values.flatten().tolist()
            core_row = core_row + core_data

            ##get fitstats
            points, vars, chi, redchi = get_fitstats(file_dir, c[0], c[1], "CORE", suffix)
            core_row = core_row + [points, vars, chi, redchi]
            data.append(core_row)
            j = j + 1

        ##flex
        if flex:
            flex_row = []
            flex = pd.read_csv(flex_file)
            flex_row.append(c[0]+'/'+c[1])
            flex_row.append("flex")
            flex_row.append("between cluster")
            flex_data = flex[flex["group"] == "all"].values.flatten().tolist()
            flex_row = flex_row + flex_data
            ##get fitstats
            points, vars, chi, redchi = get_fitstats(file_dir, c[0], c[1], "FLEX", suffix)
            flex_row = flex_row + [points, vars, chi, redchi]
            data.append(flex_row)
            j = j + 1
    print('Ran ' + str(i+j) + ' samples so far')

    if j != 0:
        between = pd.DataFrame(data,
                              columns=['cluster', 'genome', 'type', 'group','d_sample','theta_pool',
                                       'phi_pool','ratio','fbar','c','d_pool',
                                       'd_clonal','theta_s','phi_s', 'bp_analyzed', 'variables', 'chisq', 'red-chisq'])
    if i != 0 and j != 0:
        both = within.append(between)
    if i != 0 and j == 0:
        both = within
    if i == 0 and j != 0:
        both = between

    now = datetime.datetime.now()
    outpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_'+file_name+'.csv')
    both.to_csv(outpath)
    ##print a list of incomplete clusters
    repeatpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_incomplete.csv')
    repeatdf = pd.DataFrame(repeats, columns=['mate1', 'mate2','gene_type'])
    repeatdf.to_csv(repeatpath)
    ##print a list of clusters which could not be fit using lmfit (and aren't worth repeating ...)
    misfitpath = os.path.join(out_dir, now.strftime("%Y%m%d_%H%M")+'_lmfit_failed.csv')
    misfitdf = pd.DataFrame(misfits, columns=['mate1', 'mate2','gene_type'])
    misfitdf.to_csv(misfitpath)

if __name__ == "__main__":
    main()