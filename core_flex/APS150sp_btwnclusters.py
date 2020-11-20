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
job_directory = "/scratch/aps376/recombo/APS150_sp_btwnclusters"
archive = "/scratch/aps376/recombo/APS150_SP_Archive"
#scratch = os.environ['SCRATCH']
mkdir_p(archive)

"make all possible combos of clusters"

sero_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14,
             15, 17, 20, 21, 22, 23, 24, 25, 27, 28, 29, 30,
             31, 32, 34, 35, 36, 37, 38, 39, 40, 41, 43, 44, 45, 46, 47]

serocombs = combinations(sero_list, 2)
complete_combolist = []

wrkd = '/scratch/aps376/recombo/APS150_SP_Archive'
for c in serocombs:
    complete_combolist.append(('cluster'+str(c[0]), 'cluster'+str(c[1])))


gene = ['CORE', 'FLEX']

#print(str(np.arange(0,30)))
##can divide into groups of 8 for submission
count = 0
for g in gene:
    ##19 groups of 9 clusters run
    for i in np.arange(0, 83):
        combolist = complete_combolist[(i*10):((i+1)*10)]
        job_file = "%s_%s-%s.sh" % (g, str(i*10), str((i+1)*10))
        for c in combolist:
            outdir = os.path.join(archive, c[0]+'_'+c[1])
            mkdir_p(outdir)
            count = count + 1
        #os.system('cd %s' %job_directory)

        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s_%s-%s\n" % (g, str(i*10), str((i+1)*10)))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=8\n")
            fh.writelines("#SBATCH --time=12:00:00\n")
            fh.writelines("#SBATCH --mem=8GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s_%s-%s.out\n" % (g, str(i*10), str((i+1)*10)))
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
            for c in combolist:
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