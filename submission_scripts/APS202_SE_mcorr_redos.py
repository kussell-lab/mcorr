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
species = "SE_redos_v1"
date = "0521_%s" % species
jobdir = "/scratch/aps376/recombo/APS202mcorr"
archive = "/scratch/aps376/recombo/APS202_SE_Archive"
submitdir = "/scratch/aps376/recombo/APS202mcorr/%s_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS202mcorr/%s_slurm" % date
#scratch = os.environ['SCRATCH']
mkdir_p(archive)
mkdir_p(submitdir)
mkdir_p(slurmdir)

##load the csv file of failed runs ...
reruns = os.path.join(archive, "20210521_2230_unfinished_jobs.csv")
dat = pd.read_csv(reruns)
##get the bins for each
bins = np.array(dat["bin"])
##make a list of tuples for the clusters in each ...
combos = np.array(dat["cluster"])
pairs = []
for i in np.arange(0, len(combos)):
    combo = combos[i]
    if "_" in combo:
        terms = combo.split("_")
        pairs.append((terms[0], terms[1]))
    ##otherwise it's a within clade calculation ...
    else:
        pairs.append((combo, combo))

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
##
end = len(pairs)
for i in np.arange(0, end):
    b = bins[i]
    bin_archive = os.path.join(archive, "bin_" + b)
    c = pairs[i]
    count = count + 1
    #os.system('cd %s' %job_directory)
    job_file = os.path.join(submitdir, "%s_%s-%s.sh" % (b, c[0], c[1]))
    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=%s-%s_%s\n" % (c[0], c[1], b))
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=8\n")
        fh.writelines("#SBATCH --time=20:00:00\n")
        fh.writelines("#SBATCH --mem=16GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/slurm%%j_%s_%s-%s.out\n" % (slurmdir, b, c[0], c[1]))
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
        if c[0] != c[1]:
            outdir = os.path.join(archive, "bin_" + b, c[0]+'_'+c[1])
        else:
            outdir = os.path.join(archive, "bin_" + b, c[0])
        fh.writelines("cd %s\n" % outdir)
        fh.writelines("\n")
        if c[0] != c[1]:
            fh.writelines("mcorr-xmfa-2clades %s/%s/MSA_%s_%s %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT --num-boot=0 &&\n"
                          % (bin_archive, c[0], b, c[0], bin_archive, c[1], b, c[1], outdir, c[0]+'_'+c[1], b))
            fh.writelines("mcorrFitCompare %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                          % (outdir, c[0]+'_'+c[1], b, outdir, c[0]+'_'+c[1], b))
        else:
            fh.writelines("mcorr-xmfa %s/%s/MSA_%s_%s %s/%s_%s_XMFA_OUT --num-boot=0 &&\n"
                          % (bin_archive, c[0], b, c[0], outdir, c[0], b))
            fh.writelines("mcorrFitCompare %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                          % (outdir, c[0], b, outdir, c[0], b))
        fh.writelines("\n")
    os.system("sbatch %s" %job_file)
    print('submitted %s for %s-%s' %(b, c[0], c[1]))
    time.sleep(1)
print('%s total submissions' %count)