#!/usr/bin/env python
import csv
import sys
from tqdm import tqdm
import pandas as pd
from itertools import combinations
import numpy as np
from scipy import stats
# import seaborn as sns
# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt
import os

def main():

    # sero_list_file = sys.argv[1]
    # file_dir = sys.argv[2]
    # file_name = sys.argv[3]
    ##temp
    file_dir = '/Volumes/GoogleDrive/My Drive/hpc/recombo/APS148_PA_Archive/'
    file_name = 'APS148_1105_PA_all'
    out_dir = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS148_paeruginosa_analysis/'

    ##i am going to refer to these as serotypes just for convenience, but know ...
    ##they are sequence clusters!
    clusters = [1, 2, 3, 4, 5, 7, 9, 11, 12, 13, 14, 15, 16, 19, 20, 21, 22, 23, 24,
                25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 43, 44, 45,
                47, 54, 55, 56, 58, 61, 62, 63, 64, 65, 67, 68, 69, 70, 73, 74, 75, 76,
                77, 78, 79, 81, 82, 84]
    #clusters = clusters.tolist()
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

    i = 0
##for just within sero data
    for sero in tqdm(sero_list):
        core_file = file_dir + str(sero) + '/' + str(sero)+'_CORE_FIT_OUT_fit_results.csv'
        flex_file = file_dir + str(sero) + '/' + str(sero)+'_FLEX_FIT_OUT_fit_results.csv'
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
                              d, genenames, type))
        within = pd.DataFrame(all_values,
                              columns=['ST', 'phi', 'phi_median', 'phi_2.5%', 'phi_97.5%',
                                       'phi16', 'phi84',
                                       'theta', 'theta_median', 'theta_2.5%', 'theta_97.5%',
                                       'theta16', 'theta84',
                                       'd_sample', 'gene', 'type'])


    #within.to_csv(out_dir + file_name + '.csv')

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

    ##confidence level for the interval 0 to 1
    # gamma_theta = []
    # gamma_phi = []

    j = 0

    repeats = []

    for c in tqdm(combolist):
        core_file = file_dir + c[0]+'_'+c[1] + '/' + c[0]+'_'+c[1]+'_CORE_FIT_OUT_fit_results.csv'
        flex_file = file_dir + c[0]+'_'+c[1] + '/' + c[0]+'_'+c[1]+'_FLEX_FIT_OUT_fit_results.csv'
        if not os.path.exists(core_file):
            repeats.append((c[0], c[1], 'CORE'))
            continue
        elif not os.path.exists(flex_file):
            repeats.append((c[0], c[1], 'FLEX'))
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
                              d, genenames, type))
        between = pd.DataFrame(all_values,
                              columns=['ST', 'phi', 'phi_median', 'phi_2.5%', 'phi_97.5%',
                                       'phi16', 'phi84',
                                       'theta', 'theta_median', 'theta_2.5%', 'theta_97.5%',
                                       'theta16', 'theta84',
                                       'd_sample', 'gene', 'type'])
    if i != 0 and j != 0:
        both = within.append(between)
    if i != 0 and j == 0:
        both = within
    if i == 0 and j != 0:
        both = between
    both.to_csv(out_dir + file_name + '.csv')

    with open(out_dir+'repeats', 'w') as f:
        write = csv.writer(f)
        write.writerows(repeats)
   # repeats.to_csv(out_dir + 'repeats')

if __name__ == "__main__":
    main()