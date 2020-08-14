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
job_directory = "%s/APS132" %os.getcwd()
archive = "%s/Archive" %os.getcwd()
#scratch = os.environ['SCRATCH']
mkdir_p(job_directory)
mkdir_p(archive)

'make all possible combos from sero_list file, and sort them alphabetically'
sero_list_file = 'sero_list'

#read the list of serotypes
sero_list = []
with open(sero_list_file, 'r') as reader:
    for line in reader:
        sero_list.append(line.rstrip())
sero_list = sorted(sero_list)
serocombs = combinations(sero_list, 2)
combolist = []

wrkd = '/scratch/aps376/Archive'
for c in serocombs:
    combolist.append((c[0], c[1]))
 #   print(c[0]+c[1])

gene = ['CORE', 'FLEX']

for g in gene:
    for c in combolist:
        job_file = os.path.join(job_directory, "%s_%s.sh" % (g, c[0]+c[1]))
        outdir = "%s/Archive/%s_OUT" % (os.getcwd(), c[0]+c[1])
        mkdir_p(outdir)
        os.system('cd %s' %outdir)

        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s_%s\n" % (g, c[0]+c[1]))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=8\n")
            fh.writelines("#SBATCH --time=6:00:00\n")
            fh.writelines("#SBATCH --mem=32GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s/%s_%s.out\n" % (job_directory, g, c[0]+c[1]))
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("\n")
            #load modules
            fh.writelines("module purge\n")
            fh.writelines("module load git/gnu/2.16.2\n")
            fh.writelines("module load go/1.10.2\n")
            fh.writelines("module load python3/intel/3.7.3\n") ##might need to roll back to 3.6.3? unclear
            fh.writelines("module load parallel/20171022\n")
            fh.writelines("module load prokka/1.12\n")
            fh.writelines("module load muscle/intel/3.8.31\n")
            fh.writelines("module load sra-tools/2.10.5\n")
            fh.writelines("module load samtools/intel/1.6\n")
            fh.writelines("module load smalt/intel/0.7.6\n")
            fh.writelines("alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'\n")
            ##Making everything in path
            #mcorr
            fh.writelines("\n")
            fh.writelines("export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin\n")
            #ReferenceAlignmentGenerator
            fh.writelines("export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/\n")
            fh.writelines("export PATH=$PATH:~/opt/ReferenceAlignmentGenerator\n")
            fh.writelines("\n")
            fh.writelines("cd %s\n" %outdir)
            fh.writelines("\n")
            fh.writelines("mcorr-xmfa-2clades %s/%s_OUT/REFGEN_%s_%s %s/%s_OUT/REFGEN_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (wrkd, c[0], g, c[0], wrkd, c[1], g, c[1], outdir, c[0]+c[1], g))
            fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true"
                          % (outdir, c[0]+c[1], g, outdir, c[0]+c[1], g))
        os.system("sbatch %s" %job_file)
        print('submitted %s for %s' %(g, c[0]+c[1]))
        time.sleep(2)