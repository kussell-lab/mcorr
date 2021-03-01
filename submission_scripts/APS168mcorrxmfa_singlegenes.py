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
date = "0226"
jobdir = "/scratch/aps376/recombo/APS168mcorr_allseqs"
projdir="/scratch/aps376/recombo"
submitdir = "/scratch/aps376/recombo/APS168mcorr_allseqs/%s_mcorr_submissions" % date
slurmdir = "/scratch/aps376/recombo/APS168mcorr_allseqs/%s_mcorr_slurm" % date
archive = os.path.join(projdir, "APS168_SC2_Archive")
MSA= os.path.join(archive, "MSA_SC2_MASTER")
wrkd = os.path.join(archive, "genes")
#scratch = os.environ['SCRATCH']
mkdir_p(slurmdir)
mkdir_p(submitdir)

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
os.chdir(jobdir)

genes = np.arange(1, 14)
count = 0
for i in np.arange(1, 14):
    job_file = os.path.join(submitdir, "%s_mcorr.sh" % i)
    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=%s_mcorr\n" % i)
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=8\n")
        fh.writelines("#SBATCH --time=6:00:00\n")
        fh.writelines("#SBATCH --mem=30GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/slurm%%j_%s.out\n" % (slurmdir, i))
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("\n")
        #load modules
        fh.writelines("module purge\n")
        fh.writelines("module load samtools/intel/1.11\n")
        fh.writelines("module load sra-tools/2.10.9\n")
        fh.writelines("module load parallel/20201022\n")
        fh.writelines("module load python/intel/3.8.6\n")
        fh.writelines("module load smalt/intel/0.7.6\n")
        fh.writelines("module load go/1.15.7\n")
        fh.writelines("module load bowtie2/2.4.2\n")
        fh.writelines("module load bedtools/intel/2.29.2\n")
        ##Making everything in path
        #mcorr
        ##Making everything in path
        #load virtual environnment for mcorr-fit
        fh.writelines("\n")
        fh.writelines("cd /scratch/aps376/recombo\n")
        fh.writelines("source venv/bin/activate;\n")
        fh.writelines("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;\n")
        fh.writelines("export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin\n")
        fh.writelines("\n")
        fh.writelines("mcorr-xmfa %s/gene_%s %s/gene_%s_xmfa_out --show-progress --num-boot=0\n"
                      % (wrkd, i, wrkd, i))
        fh.writelines("\n")
    os.system("sbatch %s" %job_file)
    print('submitted gene %s' %i)
    count = count + 1
    time.sleep(2)
print('%s total submissions' %count)