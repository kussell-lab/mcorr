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
threshold = 95
date = "0412_%s_all" % threshold
jobdir = "/scratch/aps376/recombo/APS184mcorr"
archive = "/scratch/aps376/recombo/APS184_CJ_Archive/corethreshold%s" % threshold
submitdir = "/scratch/aps376/recombo/APS184mcorr/%s_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS184mcorr/%s_slurm" % date
#scratch = os.environ['SCRATCH']
mkdir_p(archive)
mkdir_p(submitdir)
mkdir_p(slurmdir)

"make all possible combos of clusters"

sero_list = [58, 61, 63, 73, 75, 92, 94, 95, 100, 101,
             106, 110, 112, 116, 122, 123, 125, 185, 226, 227]

serocombs = combinations(sero_list, 2)
complete_combolist = []

#wrkd = '/scratch/aps376/recombo/APS158_SP_Archive'
for c in serocombs:
    complete_combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))

for s in sero_list:
    complete_combolist.append(('cluster'+str(s), 'cluster'+str(s)))
pairs = np.array_split(complete_combolist, 150)
gene = ['CORE', 'FLEX']

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
for g in gene:
    ##will do to len(pairs) in futur
    for i in np.arange(0, len(pairs)):
        pairs_i = pairs[i]
        for c in pairs_i:
            if c[0] != c[1]:
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
            fh.writelines("#SBATCH --time=20:00:00\n")
            fh.writelines("#SBATCH --mem=16GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s/slurm%%j_%s_%s.out\n" % (slurmdir, g, i))
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("\n")
            #load modules
            fh.writelines("module purge\n")
            fh.writelines("module load go/1.15.7\n")
            fh.writelines("module load python/intel/3.8.6\n")
            ##Making everything in path
            #mcorr
            fh.writelines("\n")
            fh.writelines("export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin\n")
            #load virtual enviromment ...
            fh.writelines("cd /scratch/aps376/recombo\n")
            fh.writelines("source venv/bin/activate\n")
            fh.writelines("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;\n")
            #ReferenceAlignmentGenerator
            fh.writelines("export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/\n")
            fh.writelines("export PATH=$PATH:~/opt/ReferenceAlignmentGenerator\n")
            fh.writelines("\n")
            for c in pairs_i:
                if c[0] != c[1]:
                    outdir = os.path.join(archive, c[0]+'_'+c[1])
                else:
                    outdir = os.path.join(archive, c[0])
                fh.writelines("cd %s\n" % outdir)
                fh.writelines("\n")
                if c[0] != c[1]:
                    fh.writelines("mcorr-xmfa-2clades %s/%s/MSA_%s_%s %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                                  % (archive, c[0], g, c[0], archive, c[1], g, c[1], outdir, c[0]+'_'+c[1], g))
                    fh.writelines("mcorrFitCompare %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                                  % (outdir, c[0]+'_'+c[1], g, outdir, c[0]+'_'+c[1], g))
                else:
                    fh.writelines("mcorr-xmfa %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                                  % (archive, c[0], g, c[0], outdir, c[0], g))
                    fh.writelines("mcorrFitCompare %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                                  % (outdir, c[0], g, outdir, c[0], g))
                fh.writelines("\n")
        os.system("sbatch %s" %job_file)
        print('submitted %s job #%s' %(g, i))
        time.sleep(1)
print('%s total submissions' %count)