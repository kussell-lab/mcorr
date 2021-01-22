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
from collectresults.functions import mkdir_p, get_fitstats, get_numalns, read_lmfit_out
"""
Commandline program for collecting divergence and bootstrap interval outputs from mcorr-fit
"""

def main():
    parser = argparse.ArgumentParser(description="Collect results for many sequence clusters from just lmfit\
                                                  outputs are (1) csv file of divergences and goodness of fit stats\
                                                 (2) a list of clusters with incomplete\
                                                 jobs named 'DATETIMESTAMP_incomplete.txt")
    parser.add_argument("--file_dir", default="current", help="specify directory with output from mcorr if not current directory")
    parser.add_argument("--out_dir", default="current", help="optional out directory for output csv")
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

    if suffix == "":
        suffix = '_FIT_OUT_lmfit_report.txt'

    ##for local testing
    # file_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # out_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    # file_name = '0122_mcorr_res_test'
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

    phi = []
    theta = []
    d = []
    ##names
    genenames = []
    type = []
    seronames = []
    ##fit statistics
    alns = []
    datapoints = []
    variables = []
    chisquare = []
    reducedchisquare = []

    i = 0
    ##check if any runs didn't finish
    repeats = []
    ##ones that weren't fit by curve-fitting
    misfits = []
    ##for just within sero data
    for sero in tqdm(sero_list):
        ##check if the runs and fitting have completed
        core = True
        flex = True
        core_file = os.path.join(file_dir, str(sero), str(sero)+'_CORE_XMFA_OUT.csv')
        flex_file = os.path.join(file_dir, str(sero), str(sero)+'_FLEX_XMFA_OUT.csv')
        corestats = os.path.join(file_dir, str(sero), str(sero)+'_CORE'+suffix)
        flexstats = os.path.join(file_dir, str(sero), str(sero)+'_FLEX'+suffix)
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
            ##get phi, theta, d, goodness of fit stats, and number of alignments
            lmfitout = read_lmfit_out(file_dir, sero, "", "CORE", suffix)
            phi.append(lmfitout[6])
            theta.append(lmfitout[5])
            d.append(lmfitout[4])
            ##get fitstats
            datapoints.append(lmfitout[0])
            variables.append(lmfitout[1])
            chisquare.append(lmfitout[2])
            reducedchisquare.append(lmfitout[3])
            genenames.append('core')
            seronames.append(str(sero))
            type.append('within cluster')

            ##get avg number of alignments used to make corr profile
            avg_aln = get_numalns(file_dir, sero, "", "CORE")
            alns.append(avg_aln)

            i = i + 1

        ##flex
        if flex:
            ##get phi, theta, d, goodness of fit stats, and number of alignments
            lmfitout = read_lmfit_out(file_dir, sero, "", "FLEX", suffix)
            phi.append(lmfitout[6])
            theta.append(lmfitout[5])
            d.append(lmfitout[4])
            ##get fitstats
            datapoints.append(lmfitout[0])
            variables.append(lmfitout[1])
            chisquare.append(lmfitout[2])
            reducedchisquare.append(lmfitout[3])
            genenames.append('core')
            seronames.append(str(sero))
            type.append('within cluster')

            ##get avg number of alignments used to make corr profile
            avg_aln = get_numalns(file_dir, sero, "", "CORE")
            alns.append(avg_aln)

            i = i + 1

    if i != 0:
        all_values = list(zip(seronames, phi,
                              theta, d, alns, datapoints, variables, chisquare,
                              reducedchisquare, genenames, type))
        within = pd.DataFrame(all_values,
                              columns=['cluster', 'phi',
                                       'theta', 'd_sample', 'avg_num_alns', 'bp_analyzed', 'variables',
                                       'chi-square', 'reduced_chi-square', 'gene', 'type'])


    "for between sero data"
    phi = []
    theta = []
    d = []
    ##names
    genenames = []
    type = []
    seronames = []

    ##fit statistics
    alns = []
    datapoints = []
    variables = []
    chisquare = []
    reducedchisquare = []

    j = 0

    for c in tqdm(combolist):
        core = True
        flex = True
        core_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE_XMFA_OUT.csv')
        flex_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX_XMFA_OUT.csv')
        corestats = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE'+suffix)
        flexstats = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX'+suffix)
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
            ##get phi, theta, d, goodness of fit stats, and number of alignments
            lmfitout = read_lmfit_out(file_dir, c[0], c[1], "CORE", suffix)
            phi.append(lmfitout[6])
            theta.append(lmfitout[5])
            d.append(lmfitout[4])
            ##get fitstats
            datapoints.append(lmfitout[0])
            variables.append(lmfitout[1])
            chisquare.append(lmfitout[2])
            reducedchisquare.append(lmfitout[3])
            genenames.append('core')
            seronames.append(c[0]+'/'+c[1])
            type.append('between clusters')

            ##get avg number of alignments used to make corr profile
            avg_aln = get_numalns(file_dir, c[0], c[1], "CORE")
            alns.append(avg_aln)
            j = j + 1

        ##flex
        if flex:
            ##get phi, theta, d, goodness of fit stats, and number of alignments
            lmfitout = read_lmfit_out(file_dir, c[0], c[1], "FLEX", suffix)
            phi.append(lmfitout[6])
            theta.append(lmfitout[5])
            d.append(lmfitout[4])
            ##get fitstats
            datapoints.append(lmfitout[0])
            variables.append(lmfitout[1])
            chisquare.append(lmfitout[2])
            reducedchisquare.append(lmfitout[3])
            genenames.append('flex')
            seronames.append(c[0]+'/'+c[1])
            type.append('between clusters')

            ##get avg number of alignments used to make corr profile
            avg_aln = get_numalns(file_dir, c[0], c[1], "FLEX")
            alns.append(avg_aln)
            j = j + 1
    print('Ran ' + str(i+j) + ' samples so far')

    if j != 0:
        all_values = list(zip(seronames, phi,
                              theta, d, alns, datapoints, variables, chisquare,
                              reducedchisquare, genenames, type))
        between = pd.DataFrame(all_values,
                              columns=['cluster', 'phi',
                                       'theta', 'd_sample', 'avg_num_alns', 'bp_analyzed', 'variables',
                                       'chi-square', 'reduced_chi-square', 'gene', 'type'])
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