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
threshold = 95
date = "0318_%s_redos" % threshold
jobdir = "/scratch/aps376/recombo/APS173mcorr"
archive = "/scratch/aps376/recombo/APS173_SE_Archive/corethreshold%s" % threshold
submitdir = "/scratch/aps376/recombo/APS173mcorr/%s_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS173mcorr/%s_slurm" % date
#scratch = os.environ['SCRATCH']
mkdir_p(archive)
mkdir_p(submitdir)
mkdir_p(slurmdir)

"make a list of the redos"
lmfitfails = os.path.join(archive, "210322_1621_lmfitfailed.csv")
redos = pd.read_csv(lmfitfails)
clusterlist = np.array(redos["ID"])
genomelist = np.array(redos["genome"])

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
for i in np.arange(0, len(genomelist)):
    g = genomelist[i]
    c = clusterlist[i]
    outdir = os.path.join(archive, c)
    mkdir_p(outdir)
    count = count + 1
    job_file = os.path.join(submitdir, "%s_%s.sh" % (g, c))
    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=btwn_%s_%s\n" % (g, c))
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=8\n")
        fh.writelines("#SBATCH --time=36:00:00\n")
        fh.writelines("#SBATCH --mem=30GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/slurm%%j_%s_%s.out\n" % (slurmdir, g, c))
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
        fh.writelines("cd %s\n" % outdir)
        fh.writelines("\n")
        if len(c.split("_")) != 1:
            c = c.split("_")
            fh.writelines("mcorr-xmfa-2clades %s/%s/MSA_%s_%s %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (archive, c[0], g, c[0], archive, c[1], g, c[1], outdir, c[0]+'_'+c[1], g))
            fh.writelines("mcorrFitOne %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                          % (outdir, c[0]+'_'+c[1], g, outdir, c[0]+'_'+c[1], g))
        else:
            c = c.split("_")
            fh.writelines("mcorr-xmfa %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (archive, c[0], g, c[0], outdir, c[0], g))
            fh.writelines("mcorrFitOne %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                          % (outdir, c[0], g, outdir, c[0], g))
        fh.writelines("\n")
    os.system("sbatch %s" %job_file)
    print('submitted %s job #%s' %(g, i))
    time.sleep(1)
print('%s total submissions' %count)