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
from collectresults.functions import mkdir_p, get_fitstats
"""
Commandline program for collecting divergence and bootstrap interval outputs from mcorr-fit
"""

# def mkdir_p(dir):
#     'make a directory if doesnt exist'
#     if not os.path.exists(dir):
#         os.mkdir(dir)

def main():
    parser = argparse.ArgumentParser(description="Collect results for many sequence clusters from mcorr-fit\n"+
                                     "outputs are (1) csv file of divergences, (2) a list of clusters with incomplete\n"+
                                                 " jobs named 'DATETIMESTAMP_incomplete.txt'")
    parser.add_argument("--file_dir", default="current", help="specify directory with output from mcorr if not current directory")
    parser.add_argument("--out_dir", default="current", help="optional output directory for divergence csv")
    parser.add_argument("stats_sheet", help=" name of stats csv file from clusterSequences.py which has the list of clusters")
    parser.add_argument("file_name", help="prefix for resultant csv file (timestamp will be included)")

    ##define commandline args as variables
    args = parser.parse_args()
    file_dir = args.file_dir
    file_name = args.file_name
    stats_sheet = args.stats_sheet
    out_dir = args.out_dir

    ##for local testing
    #file_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    #out_dir = '/Volumes/aps_timemachine/recombo/APS160.5_lmfit'
    #file_name = 'APS160.5_mcorr_res_test'
    #clusters = [8, 221]
    # stats_sheet = "20th_percentile_>=10strains_stats.csv"

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
    ##bootstrap medians
    phib = []
    thetab = []
    ##bootstrap CIs
    phi2pt5 = []
    phi97pt5 = []
    theta2pt5 = []
    theta97pt5 = []
    phi16 = []
    phi84 = []
    theta16 = []
    theta84 = []
    ##names
    genenames = []
    type = []
    seronames = []
    ##fit statistics
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
        missing = False
        core_file = os.path.join(file_dir, str(sero), str(sero)+'_CORE_FIT_OUT_fit_results.csv')
        flex_file = os.path.join(file_dir, str(sero), str(sero)+'_FLEX_FIT_OUT_fit_results.csv')
        corestats = os.path.join(file_dir, str(sero), str(sero)+'_CORE_FIT_OUT_lmfit_report.txt')
        flexstats = os.path.join(file_dir, str(sero), str(sero)+'_FLEX_FIT_OUT_lmfit_report.txt')
        if not os.path.exists(core_file):
            repeats.append((sero, sero, 'CORE'))
            missing = True
        if not os.path.exists(flex_file):
            repeats.append((sero, sero, 'FLEX'))
            missing = True
        if not os.path.exists(corestats):
            misfits.append((sero, sero, 'CORE'))
            missing = True
        if not os.path.exists(flexstats):
            misfits.append((sero, sero, 'FLEX'))
            missing = True
        if missing:
            continue
        core = pd.read_csv(core_file)
        flex = pd.read_csv(flex_file)

        ##core
        phi.append(core['phi_pool'][0])
        theta.append(core['theta_pool'][0])
        d.append(core['d_sample'][0])
        genenames.append('core')
        seronames.append(str(sero))
        type.append('within serotype')

        ##get bootstrap medians and CIs
        phi_pool = core['phi_pool'][1:]
        theta_pool = core['theta_pool'][1:]
        phib.append(np.median(phi_pool))
        thetab.append(np.median(theta_pool))
        phi2pt5.append(np.percentile(phi_pool, 2.5))
        phi97pt5.append(np.percentile(phi_pool, 97.5))
        theta2pt5.append(np.percentile(theta_pool, 2.5))
        theta97pt5.append(np.percentile(theta_pool, 97.5))
        phi16.append(np.percentile(phi_pool, 16))
        phi84.append(np.percentile(phi_pool, 84))
        theta16.append(np.percentile(theta_pool, 16))
        theta84.append(np.percentile(theta_pool, 84))

        ##get fitstats
        points, vars, chi, redchi = get_fitstats(file_dir, sero, "", "CORE")
        datapoints.append(points)
        variables.append(vars)
        chisquare.append(chi)
        reducedchisquare.append(redchi)

        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('within serotype')
        seronames.append(str(sero))

        ##get bootstrap medians
        phi_pool = flex['phi_pool'][1:]
        theta_pool = flex['theta_pool'][1:]
        phib.append(np.median(phi_pool))
        thetab.append(np.median(theta_pool))
        phi2pt5.append(np.percentile(phi_pool, 2.5))
        phi97pt5.append(np.percentile(phi_pool, 97.5))
        theta2pt5.append(np.percentile(theta_pool, 2.5))
        theta97pt5.append(np.percentile(theta_pool, 97.5))
        phi16.append(np.percentile(phi_pool, 16))
        phi84.append(np.percentile(phi_pool, 84))
        theta16.append(np.percentile(theta_pool, 16))
        theta84.append(np.percentile(theta_pool, 84))

        ##get fitstats
        points, vars, chi, redchi = get_fitstats(file_dir, sero, "", "FLEX")
        datapoints.append(points)
        variables.append(vars)
        chisquare.append(chi)
        reducedchisquare.append(redchi)


        ##count how many you've done so far
        if os.path.exists(core_file) and os.path.exists(flex_file):
            i = i + 2
        else:
            i = i + 1
    if i != 0:
        all_values = list(zip(seronames, phi, phib, phi2pt5, phi97pt5,
                              phi16, phi84,
                              theta, thetab, theta2pt5, theta97pt5,
                              theta16, theta84,
                              d, datapoints, variables, chisquare, reducedchisquare, genenames, type))
        within = pd.DataFrame(all_values,
                              columns=['cluster', 'phi', 'phi_median', 'phi_2.5%', 'phi_97.5%',
                                       'phi16', 'phi84',
                                       'theta', 'theta_median', 'theta_2.5%', 'theta_97.5%',
                                       'theta16', 'theta84', 'd_sample', 'bp analyzed', 'variables',
                                       'chi-square', 'reduced chi-square', 'gene', 'type'])


    "for between sero data"
    phi = []
    theta = []
    d = []
    ##names
    genenames = []
    type = []
    seronames = []

    ##bootstrap medians
    phib = []
    thetab = []
    ##bootstrap CIs
    phi2pt5 = []
    phi97pt5 = []
    theta2pt5 = []
    theta97pt5 = []
    phi16 = []
    phi84 = []
    theta16 = []
    theta84 = []

    ##fit statistics
    datapoints = []
    variables = []
    chisquare = []
    reducedchisquare = []

    j = 0

    for c in tqdm(combolist):
        missing = False
        core_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE_FIT_OUT_fit_results.csv')
        flex_file = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX_FIT_OUT_fit_results.csv')
        corestats = os.path.join(file_dir, c[0]+'_'+c[1], c[0]+'_'+c[1]+'_CORE_FIT_OUT_lmfit_report.txt')
        flexstats = os.path.join(file_dir,  c[0]+'_'+c[1], c[0]+'_'+c[1]+'_FLEX_FIT_OUT_lmfit_report.txt')
        if not os.path.exists(core_file):
            repeats.append((c[0], c[1], 'CORE'))
            missing = True
        if not os.path.exists(flex_file):
            repeats.append((c[0], c[1], 'FLEX'))
            missing = True
        if not os.path.exists(corestats):
            misfits.append((c[0], c[1], 'CORE'))
            missing = True
        if not os.path.exists(flexstats):
            misfits.append((c[0], c[1], 'FLEX'))
            missing = True
        if missing:
            continue
        core = pd.read_csv(core_file)
        flex = pd.read_csv(flex_file)
        ##core
        phi.append(core['phi_pool'][0])
        theta.append(core['theta_pool'][0])
        d.append(core['d_sample'][0])
        genenames.append('core')
        type.append('between serotype')
        seronames.append(c[0]+'/'+c[1])

        ##get bootstrap medians
        phi_pool = core['phi_pool'][1:]
        theta_pool = core['theta_pool'][1:]
        phib.append(np.median(phi_pool))
        thetab.append(np.median(theta_pool))
        phi2pt5.append(np.percentile(phi_pool, 2.5))
        phi97pt5.append(np.percentile(phi_pool, 97.5))
        theta2pt5.append(np.percentile(theta_pool, 2.5))
        theta97pt5.append(np.percentile(theta_pool, 97.5))
        phi16.append(np.percentile(phi_pool, 16))
        phi84.append(np.percentile(phi_pool, 84))
        theta16.append(np.percentile(theta_pool, 16))
        theta84.append(np.percentile(theta_pool, 84))

        ##get goodness of fit statistics
        points, vars, chi, redchi = get_fitstats(file_dir, c[0], c[1], "CORE")
        datapoints.append(points)
        variables.append(vars)
        chisquare.append(chi)
        reducedchisquare.append(redchi)
        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('between serotype')
        seronames.append(c[0]+'/'+c[1])

        ##get bootstrap medians
        phi_pool = flex['phi_pool'][1:]
        theta_pool = flex['theta_pool'][1:]
        phib.append(np.median(phi_pool))
        thetab.append(np.median(theta_pool))
        phi2pt5.append(np.percentile(phi_pool, 2.5))
        phi97pt5.append(np.percentile(phi_pool, 97.5))
        theta2pt5.append(np.percentile(theta_pool, 2.5))
        theta97pt5.append(np.percentile(theta_pool, 97.5))
        phi16.append(np.percentile(phi_pool, 16))
        phi84.append(np.percentile(phi_pool, 84))
        theta16.append(np.percentile(theta_pool, 16))
        theta84.append(np.percentile(theta_pool, 84))

        ##get goodness of fit statistics
        points, vars, chi, redchi = get_fitstats(file_dir, c[0], c[1], "FLEX")
        datapoints.append(points)
        variables.append(vars)
        chisquare.append(chi)
        reducedchisquare.append(redchi)

        if os.path.exists(core_file) and os.path.exists(flex_file):
            j = j + 2
        else:
            j = j + 1
    print('Ran ' + str(i+j) + ' samples so far')
    ##make sure you're not appending

    if j != 0:
        all_values = list(zip(seronames, phi, phib, phi2pt5, phi97pt5,
                              phi16, phi84,
                              theta, thetab, theta2pt5, theta97pt5,
                              theta16, theta84,
                              d, datapoints, variables, chisquare, reducedchisquare, genenames, type))
        between = pd.DataFrame(all_values,
                              columns=['cluster', 'phi', 'phi_median', 'phi_2.5%', 'phi_97.5%',
                                       'phi16', 'phi84',
                                       'theta', 'theta_median', 'theta_2.5%', 'theta_97.5%',
                                       'theta16', 'theta84', 'd_sample', 'bp analyzed', 'variables',
                                       'chi-square', 'reduced chi-square', 'gene', 'type'])
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