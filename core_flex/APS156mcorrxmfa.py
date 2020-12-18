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
jobdir = "/scratch/aps376/recombo/APS156mcorr-xmfa/1207mcorr"
archive = "/scratch/aps376/recombo/APS156_SP_Archive"
outdir = "/scratch/aps376/recombo/APS156_SP_Archive/mcorr_out"
wrkd = "/scratch/aps376/recombo/APS156_SP_Archive/genes"
#scratch = os.environ['SCRATCH']
mkdir_p(outdir)
mkdir_p(jobdir)

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
os.chdir(jobdir)

genes = np.arange(1, 2020)
gene_groups = np.array_split(genes, 100)
count = 0
for i in np.arange(0,len(gene_groups)):
        job_file = "group_%s.sh" % i
        genes_i = gene_groups[i]
        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=mcorr_%s\n" % i)
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=4\n")
            fh.writelines("#SBATCH --time=6:00:00\n")
            fh.writelines("#SBATCH --mem=4GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=mcorr_%s_slurm%%j.out\n" % i)
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
            for gene in genes_i:
                count = count + 1
                fh.writelines("\n")
                fh.writelines("mcorr-xmfa %s/gene_%s %s/gene_%s_xmfa_out --show-progress --num-boot=0\n"
                              % (wrkd, gene, outdir, gene))
                fh.writelines("\n")
        os.system("sbatch %s" %job_file)
        print('submitted group %s' %i)
        time.sleep(2)
print('%s total submissions' %count)