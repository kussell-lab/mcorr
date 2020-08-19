#!/usr/bin/env python
import sys
from tqdm import tqdm
import pandas as pd
from itertools import combinations
# import seaborn as sns
# import matplotlib as mpl
# mpl.use('Agg')
# import matplotlib.pyplot as plt
import os

def main():

    sero_list_file = sys.argv[1]
    file_dir = sys.argv[2]
    ##temp
    #file_dir = '/Volumes/GoogleDrive/My Drive/200818_Archive/'

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
        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('within serotype')
        seronames.append(str(sero))
        if os.path.exists(core_file) and os.path.exists(flex_file):
            i = i + 2
        else:
            i = i + 1
    if i != 0:
        all_values = list(zip(seronames, phi, theta, d, genenames, type))
        within = pd.DataFrame(all_values, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])
    #within.to_csv(file_dir + 'corevflex.csv') #, sep='\t')
####for between sero data
    phi = []
    theta = []
    d = []
    ##names
    genenames = []
    type = []
    seronames = []
    j = 0
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
        ##flex
        phi.append(flex['phi_pool'][0])
        theta.append(flex['theta_pool'][0])
        d.append(flex['d_sample'][0])
        genenames.append('flex')
        type.append('between serotype')
        seronames.append(c[0]+c[1])
        if os.path.exists(core_file) and os.path.exists(flex_file):
            j = j + 2
        else:
            j = j + 1
    print('Ran ' + str(i+j) + ' samples so far')

    if j != 0:
        all_values = list(zip(seronames, phi, theta, d, genenames, type))
        between = pd.DataFrame(all_values, columns=['sero', 'phi', 'theta', 'd_sample', 'gene', 'type'])
    if i != 0 and j != 0:
        both = within.append(between)
    if i != 0 and j == 0:
        both = within
    if i == 0 and j != 0:
        both = between
    both.to_csv(file_dir + 'corevflex_w_diversity.csv') #, sep='\t')
    ####
    #plotting currently not working on HPC (which is fine), will ask them to update
    #python 3.7 with the necessary packages
    # supescool=['#d9544d','#3b5b92']
    # sns.set_style("darkgrid")
    # sns.set(font_scale = 1.5)
    # fig_dims = (7, 6)
    # fig, ax = plt.subplots(figsize=fig_dims)
    # sns.scatterplot(x='theta', y = 'phi', hue='gene', style='type', data=both,
    #                 s = 50,
    #                 alpha = 0.7, edgecolor = 'k', palette = supescool, ax = ax)
    #
    # ax.set_xlabel(r'$\theta$')
    # ax.set_ylabel('$\phi$')
    # plt.savefig(file_dir + 'corevflex.pdf', bbox_inches="tight")
if __name__ == "__main__":
    main()