#!/usr/bin/env python
import os
import time
from itertools import combinations
import numpy as np
import pandas as pd


def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)

#define directories
date = "0109"
jobdir = "/scratch/aps376/recombo/APS158_sp_btwnclusters"
archive = "/scratch/aps376/recombo/APS158_SP_Archive"
submitdir = "/scratch/aps376/recombo/APS158_sp_btwnclusters/%s_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS158_sp_btwnclusters/%s_slurm" % date
clusterlist = "/scratch/aps376/recombo/APS158_SP_Archive/20210109_1208_incomplete.csv"
#scratch = os.environ['SCRATCH']
mkdir_p(archive)
mkdir_p(submitdir)
mkdir_p(slurmdir)

#clusterlist = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS150_SP_analysis/20210109_1204_incomplete.csv'
"load cluster list a list of tuples"
clusterdf = pd.read_csv(clusterlist, index_col=0)
records = clusterdf.to_records(index=False)
clusters = list(records)
#print(clusters)

wrkd = '/scratch/aps376/recombo/APS158_SP_Archive'

# sero_list = [8, 9, 27, 75, 83, 85, 89, 94, 99, 106, 110, 111, 112, 133,
#              136, 140, 141, 142, 149, 152, 155, 158, 159, 161, 162, 163,
#              165, 169, 170, 171, 173, 174, 178, 180, 183, 191, 196, 198,
#              199, 213, 216, 218, 220, 221, 228, 229]
#
# serocombs = combinations(sero_list, 2)
# complete_combolist = []
#
# wrkd = '/scratch/aps376/recombo/APS158_SP_Archive'
# for c in serocombs:
#     complete_combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))
#
# pairs = np.array_split(complete_combolist, 250)
# gene = ['CORE', 'FLEX']

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
for c in clusters:
    ##will do to len(pairs) in futur
    outdir = os.path.join(archive, c[0]+'_'+c[1])
    mkdir_p(outdir)
    count = count + 1
    #os.system('cd %s' %job_directory)
    job_file = os.path.join(submitdir, "%s-%s_%s.sh" % (c[0], c[1], c[2]))
    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=%s-%s_%s\n" % (c[0], c[1], c[2]))
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=4\n")
        fh.writelines("#SBATCH --time=12:00:00\n")
        fh.writelines("#SBATCH --mem=4GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/slurm%%j_%s-%s_%s.out\n" % (slurmdir, c[0], c[1], c[2]))
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("\n")
        #load modules
        fh.writelines("module purge\n")
        fh.writelines("module load go/1.15.2\n")
        fh.writelines("module load python/intel/3.8.6\n")
        ##Making everything in path
        #mcorr
        fh.writelines("\n")
        fh.writelines("export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin\n")
        #ReferenceAlignmentGenerator
        fh.writelines("export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/\n")
        fh.writelines("export PATH=$PATH:~/opt/ReferenceAlignmentGenerator\n")
        fh.writelines("\n")
        outdir = os.path.join(archive, c[0]+'_'+c[1])
        fh.writelines("cd %s\n" %outdir)
        fh.writelines("\n")
        fh.writelines("mcorr-xmfa-2clades %s/%s/MSA_%s_%s %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (wrkd, c[0], c[2], c[0], wrkd, c[1], c[2], c[1], outdir, c[0]+'_'+c[1], c[2]))
        fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                      % (outdir, c[0]+'_'+c[1], c[2], outdir, c[0]+'_'+c[1], c[2]))
    os.system("sbatch %s" %job_file)
    print('submitted %s-%s for %s' %(c[0], c[1], c[2]))
    time.sleep(2)
print('%s total submissions' %count)