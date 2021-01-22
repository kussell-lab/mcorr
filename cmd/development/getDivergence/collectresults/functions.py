#!/usr/bin/env python3
import os
import tempfile
import pandas as pd
import numpy as np

def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)
def get_fitstats(file_dir, name1, name2, gene, suffix):
    'get goodness of fit stats'
    if suffix == "":
        suffix = '_FIT_OUT_lmfit_report.txt'
    if name2 == "":
        statsfile = os.path.join(file_dir, str(name1), str(name1)+'_'+gene+suffix)
    else:
        statsfile = os.path.join(file_dir, name1+'_'+name2, name1+'_'+name2+'_'+gene+suffix)
    stats = open(statsfile)
    for i, line in enumerate(stats):
        if i == 3:
            terms = line.rstrip().split(" ")
            datapoints = terms[len(terms)-1]
        if i == 4:
            terms = line.rstrip().split(" ")
            variables = terms[len(terms)-1]
        if i == 5:
            terms = line.rstrip().split(" ")
            chisquare = terms[len(terms)-1]
        if i == 6:
            terms = line.rstrip().split(" ")
            reducedchisquare = terms[len(terms)-1]
        if i == 7:
            break
    return datapoints, variables, chisquare, reducedchisquare

def get_numalns(file_dir, name1, name2, gene):
    'get the average number of alignments used to calculate the correlation profile'
    if name2 == "":
        corrfile = os.path.join(file_dir, str(name1), str(name1)+'_'+gene+'_XMFA_OUT.csv')
    else:
        corrfile = os.path.join(file_dir, name1+'_'+name2, name1+'_'+name2+'_'+gene+'_XMFA_OUT.csv')
    corr = pd.read_csv(corrfile, skiprows=6)
    allcorr = corr[corr['b'] == "all"]
    allcorr = allcorr[allcorr['t'] == 'P2']
    numaln = allcorr['n'].to_numpy()
    avg_aln = np.mean(numaln)
    return avg_aln

def read_lmfit_out(file_dir, name1, name2, gene, suffix):
    'read lmfit output and return phi_pool, theta_pool, d_sample and goodness of fit stats'
    if suffix == "":
        suffix = '_FIT_OUT_lmfit_report.txt'
    if name2 == "":
        statsfile = os.path.join(file_dir, str(name1), str(name1)+'_'+gene+suffix)
    else:
        statsfile = os.path.join(file_dir, name1+'_'+name2, name1+'_'+name2+'_'+gene+suffix)
    stats = open(statsfile)
    for i, line in enumerate(stats):
        terms = line.strip().split(" ")
        if len(terms) > 1:
            if terms[1] == "data":
                datapoints = terms[len(terms)-1]
            elif terms[1] == "variables":
                variables = terms[len(terms)-1]
            elif terms[0] == "chi-square":
                chisquare = terms[len(terms)-1]
            elif terms[0] == "reduced":
                reducedchisquare = terms[len(terms)-1]
            elif terms[0] == "ds:":
                d_sample = terms[len(terms)-2]
            elif terms[0] == "thetaS:":
                theta_sample = terms[2]
            elif terms[0] == "f:":
                f = terms[7]
            elif terms[0] == "phiS:":
                phi_sample = terms[4]
            elif terms[0] == "thetaP:":
                theta_pool = terms[2]
            elif terms[0] == "phiP:":
                phi_pool = terms[4]
                break

    return datapoints, variables, chisquare, reducedchisquare, \
           d_sample, theta_sample, f, phi_sample, theta_pool, phi_pool
