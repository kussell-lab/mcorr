#!/usr/bin/env python
import os
import time
from itertools import combinations
import numpy as np


def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)

#define directories
date = "0113"
jobdir = "/scratch/aps376/recombo/APS159_ng_btwnclusters"
archive = "/scratch/aps376/recombo/APS159_NG_Archive"
submitdir = "/scratch/aps376/recombo/APS159_ng_btwnclusters/%s_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS159_ng_btwnclusters/%s_slurm" % date
#scratch = os.environ['SCRATCH']
mkdir_p(archive)
mkdir_p(submitdir)
mkdir_p(slurmdir)

"make all possible combos of clusters"

sero_list = [3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15, 16, 17, 19, 20, 21, 22,
             23, 26, 28, 29, 31, 32, 33, 34, 36, 37, 39, 40, 41, 43, 46, 48,
             49, 50, 51, 52, 53, 54, 57, 58, 59, 60, 61, 62]

serocombs = combinations(sero_list, 2)
complete_combolist = []

wrkd = '/scratch/aps376/recombo/APS159_NG_Archive'
for c in serocombs:
    complete_combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))

pairs = np.array_split(complete_combolist, 100)
gene = ['CORE', 'FLEX']

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
for g in gene:
    ##will do to len(pairs) in futur
    for i in np.arange(0, len(pairs)):
        pairs_i = pairs[i]
        # combolist = complete_combolist[(i*10):((i+1)*10)]
        # job_file = "%s_%s-%s.sh" % (g, str(i*10), str((i+1)*10))
        for c in pairs_i:
            outdir = os.path.join(archive, c[0]+'_'+c[1])
            mkdir_p(outdir)
            count = count + 1
        #os.system('cd %s' %job_directory)
        job_file = os.path.join(submitdir, "%s_group_%s.sh" % (g, i))
        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=btwn_%s_%s\n" % (g, str(i)))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=4\n")
            fh.writelines("#SBATCH --time=12:00:00\n")
            fh.writelines("#SBATCH --mem=4GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s/slurm%%j_%s_%s.out\n" % (slurmdir, g, i))
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
            for c in pairs_i:
                outdir = os.path.join(archive, c[0]+'_'+c[1])
                fh.writelines("cd %s\n" %outdir)
                fh.writelines("\n")
                fh.writelines("mcorr-xmfa-2clades %s/%s/MSA_%s_%s %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                              % (wrkd, c[0], g, c[0], wrkd, c[1], g, c[1], outdir, c[0]+'_'+c[1], g))
                fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                              % (outdir, c[0]+'_'+c[1], g, outdir, c[0]+'_'+c[1], g))
                fh.writelines("\n")
        os.system("sbatch %s" %job_file)
        print('submitted %s job #%s' %(g, i))
        time.sleep(2)
print('%s total submissions' %count)