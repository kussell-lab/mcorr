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
