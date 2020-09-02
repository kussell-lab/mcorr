#!/usr/bin/env python
import sys
from tqdm import tqdm
import pandas as pd
from itertools import combinations
import numpy as np
# import seaborn as sns
# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt
import os

def main():

    sero_list_file = sys.argv[1]
    file_dir = sys.argv[2]
    file_name = sys.argv[3]
    ##temp
    file_dir = '/Volumes/GoogleDrive/My Drive/hpc/APS135/APS135_200823_Archive/'

    #read the list of serotypes
    sero_list = []
    with open(sero_list_file, 'r') as reader:
        for line in reader:
            sero_list.append(line.rstrip())
    sero_list = sorted(sero_list)
    serocombs = combinations(sero_list, 2)
    combolist = []
    for c in serocombs:
        combolist.append((c[0], c[1]))

    phi = []
    theta = []
    d = []

    ##names
    genenames = []
    type = []
    seronames = []

    ##filtered by 95% bootstrap CIs
    phiF = []
    thetaF = []
    dF = []
    genenamesF = []
    typeF = []
    seronamesF = []
    f = 0
    i = 0
##for just within sero data
    for sero in tqdm(sero_list):
        core_file = file_dir + str(sero)+'_OUT/'+str(sero)+'_CORE_FIT_OUT_fit_results.csv'
        flex_file = file_dir + str(sero)+'_OUT/'+str(sero)+'_FLEX_FIT_OUT_fit_results.csv'
        if not os.path.exists(core_file):
            continue
        elif not os.path.exists(flex_file):
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

        ##filter core values by bootstrap CI
        phi_pool = core['phi_pool'][1:]
        theta_pool = core['theta_pool'][1:]
        phic = core['phi_pool'][0]
        thetac = core['theta_pool'][0]
        if np.percentile(phi_pool, 2.5) < phic and phic < np.percentile(phi_pool, 97.5):
            if np.percentile(theta_pool, 2.5) < thetac and thetac < np.percentile(theta_pool, 97.5):
                phiF.append(core['phi_pool'][0])
                thetaF.append(core['theta_pool'][0])
                dF.append(core['d_sample'][0])
                genenamesF.append('core')
                seronamesF.append(str(sero))
                typeF.append('within serotype')
                f = f + 1
        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('within serotype')
        seronames.append(str(sero))

        ##filter flex values by bootstrap CI
        phi_pool = flex['phi_pool'][1:]
        theta_pool = flex['theta_pool'][1:]
        phifl = flex['phi_pool'][0]
        thetafl = flex['theta_pool'][0]
        if np.percentile(phi_pool, 2.5) < phifl and phifl < np.percentile(phi_pool, 97.5):
            if np.percentile(theta_pool, 2.5) < thetafl and thetafl < np.percentile(theta_pool, 97.5):
                phiF.append(flex['phi_pool'][0])
                thetaF.append(flex['theta_pool'][0])
                dF.append(flex['d_sample'][0])
                genenamesF.append('flex')
                seronamesF.append(str(sero))
                typeF.append('within serotype')
                f = f + 1

        ##count how many you've done so far
        if os.path.exists(core_file) and os.path.exists(flex_file):
            i = i + 2
        else:
            i = i + 1
    if i != 0:
        all_values = list(zip(seronames, phi, theta, d, genenames, type))
        within = pd.DataFrame(all_values, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])
        ### now data which has been filtered by 95% CIs
    if f != 0:
        filtered = list(zip(seronamesF, phiF, thetaF, dF, genenamesF, typeF))
        withinF = pd.DataFrame(filtered, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])

    "for between sero data"
    phi = []
    theta = []
    d = []
    ##names
    genenames = []
    type = []
    seronames = []

    ##filtered by 95% bootstrap CIs
    phiF = []
    thetaF = []
    dF = []
    genenamesF = []
    typeF = []
    seronamesF = []

    j = 0
    f1 = 0
    for c in combolist:
        core_file = file_dir + c[0]+c[1]+'_OUT/'+c[0]+c[1]+'_CORE_FIT_OUT_fit_results.csv'
        flex_file = file_dir + c[0]+c[1]+'_OUT/'+c[0]+c[1]+'_FLEX_FIT_OUT_fit_results.csv'
        if not os.path.exists(core_file):
            continue
        elif not os.path.exists(flex_file):
            continue
        core = pd.read_csv(core_file)
        flex = pd.read_csv(flex_file)
        ##core
        phi.append(core['phi_pool'][0])
        theta.append(core['theta_pool'][0])
        d.append(core['d_sample'][0])
        genenames.append('core')
        type.append('between serotype')
        seronames.append(c[0]+c[1])

        ##filter core values by bootstrap CI
        phi_pool = core['phi_pool'][1:]
        theta_pool = core['theta_pool'][1:]
        phic = core['phi_pool'][0]
        thetac = core['theta_pool'][0]
        if np.percentile(phi_pool, 2.5) < phic and phic < np.percentile(phi_pool, 97.5):
            if np.percentile(theta_pool, 2.5) < thetac and thetac < np.percentile(theta_pool, 97.5):
                phiF.append(core['phi_pool'][0])
                thetaF.append(core['theta_pool'][0])
                dF.append(core['d_sample'][0])
                genenamesF.append('core')
                seronamesF.append(str(sero))
                typeF.append('within serotype')
                f1 = f1 + 1

        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('between serotype')
        seronames.append(c[0]+c[1])

        ##filter flex values by bootstrap CI
        phi_pool = flex['phi_pool'][1:]
        theta_pool = flex['theta_pool'][1:]
        phifl = flex['phi_pool'][0]
        thetafl = flex['theta_pool'][0]
        if np.percentile(phi_pool, 2.5) < phifl and phifl < np.percentile(phi_pool, 97.5):
            if np.percentile(theta_pool, 2.5) < thetafl and thetafl < np.percentile(theta_pool, 97.5):
                phiF.append(flex['phi_pool'][0])
                thetaF.append(flex['theta_pool'][0])
                dF.append(flex['d_sample'][0])
                genenamesF.append('flex')
                seronamesF.append(str(sero))
                typeF.append('within serotype')
                f1 = f1 + 1

        if os.path.exists(core_file) and os.path.exists(flex_file):
            j = j + 2
        else:
            j = j + 1
    print('Ran ' + str(i+j) + ' samples so far')
    print('Ran ' + str(f+f1) + ' filtered samples so far')
    ##make sure you're not appending

    if j != 0:
        all_values = list(zip(seronames, phi, theta, d, genenames, type))
        between = pd.DataFrame(all_values, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])
    if i != 0 and j != 0:
        both = within.append(between)
    if i != 0 and j == 0:
        both = within
    if i == 0 and j != 0:
        both = between
    both.to_csv(file_dir + file_name + '.csv') #, sep='\t')
    ####for the 95% bootstrap CI filtered data ...
    if f1 != 0:
        filtered1 = list(zip(seronamesF, phiF, thetaF, dF, genenamesF, typeF))
        betweenF = pd.DataFrame(filtered1, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])
    if f != 0 and f1 != 0:
        bothF = withinF.append(betweenF)
    if f != 0 and f1 == 0:
        both = within
    if f == 0 and f1 != 0:
        both = between
    bothF.to_csv(file_dir + file_name + '_filtered.csv') #, sep='\t')

if __name__ == "__main__":
    main()