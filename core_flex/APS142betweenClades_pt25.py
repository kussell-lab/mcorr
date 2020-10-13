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
job_directory = "%s/APS142_pt25_btwnclades" %os.getcwd()
archive = "%s/APS142_ngs_pt25_Archive" %os.getcwd()
msadir = "%s/CF_cutoff0.25" %os.getcwd()
#scratch = os.environ['SCRATCH']
mkdir_p(job_directory)
mkdir_p(archive)

'make all possible combos of clusters, and sort them alphanumerically'
#sero_list_file = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS137_ngs/ST_list'

##cluster names

clusters = [1, 4, 5, 6, 7, 9, 10, 11, 12, 13, 14, 15, 16, 20, 21, 22, 23, 24, 25,
            26, 27, 30, 32, 33, 34, 35, 36, 37, 38, 39, 43, 44, 45]
serocombs = combinations(clusters, 2)
combolist = []

wrkd = '/scratch/aps376/APS142_ngs_pt25_Archive'
for c in serocombs:
    combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))

gene = ['CORE', 'FLEX']

#print(len(combolist))
##can divide into groups of 8 for submission
count = 0

for c in combolist:
        job_file = os.path.join(job_directory, "%s-%s.sh" % (c[0], c[1]))
        outdir = archive
        os.system('cd %s' %job_directory)

        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s-%s\n" % (c[0], c[1]))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=4\n")
            fh.writelines("#SBATCH --time=4:00:00\n")
            fh.writelines("#SBATCH --mem=8GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s/%s-%s.out\n" % (job_directory, c[0], c[1]))
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
            fh.writelines("mcorr-xmfa-2clades %s/MSA_%s_%s %s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                              % (msadir, gene[0], c[0], msadir, gene[0], c[1], outdir, c[0]+'_'+c[1], gene[0]))
            fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                              % (outdir, c[0]+'_'+c[1], gene[0], outdir, c[0]+'_'+c[1], gene[0]))
            fh.writelines("mcorr-xmfa-2clades %s/MSA_%s_%s %s/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (msadir, gene[1], c[0], msadir, gene[1], c[1], outdir, c[0]+'_'+c[1], gene[1]))
            fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                          % (outdir, c[0]+'_'+c[1], gene[1], outdir, c[0]+'_'+c[1], gene[1]))
        os.system("sbatch %s" %job_file)
        count = count + 1
        print('submitted %s and %s for %s and %s' %(gene[0], gene[1], c[0], c[1]))
        time.sleep(2)
print('%s total submissions' %str(count))