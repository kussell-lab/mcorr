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
job_directory = "%s/APS137_btwnclades" %os.getcwd()
archive = "%s/APS137_Archive" %os.getcwd()
#scratch = os.environ['SCRATCH']
mkdir_p(job_directory)
mkdir_p(archive)

'make all possible combos from sero_list file, and sort them alphabetically'
sero_list_file = '/home/aps376/APS137_ngs/ST_list'
#sero_list_file = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS137_ngs/ST_list'

#read the list of serotypes
sero_list = []
with open(sero_list_file, 'r') as reader:
    for line in reader:
        sero_list.append(line.rstrip())
sero_list = sorted(sero_list)
serocombs = combinations(sero_list, 2)
complete_combolist = []

wrkd = '/scratch/aps376/APS137_Archive'
for c in serocombs:
    complete_combolist.append((str(c[0]), str(c[1])))
 #   print(c[0]+c[1])

#combolist = [(str(11172), str(1588))]

gene = ['CORE', 'FLEX']

#print(len(combolist))
##can divide into groups of 8 for submission
count = 0
for g in gene:
    ##will go to 66 on cluster
    for i in np.arange(0, 66):
        combolist = complete_combolist[(i*8):((i+1)*8)]
        job_file = os.path.join(job_directory, "%s_ST-%s-%s.sh" % (g, str(i*8), str((i+1)*8)))
        for c in combolist:
            outdir = "%s/APS137_Archive/%s_OUT" % (os.getcwd(), 'ST-'+c[0]+'_ST-'+c[1])
            mkdir_p(outdir)
            count = count + 1
        os.system('cd %s' %job_directory)

        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s_ST-%s-%s\n" % (g, str(i*8), str((i+1)*8)))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=8\n")
            fh.writelines("#SBATCH --time=6:00:00\n")
            fh.writelines("#SBATCH --mem=8GB\n")
            fh.writelines("#SBATCH --mail-type=END,FAIL\n")
            fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
            fh.writelines("#SBATCH --output=%s/%s_ST-%s-%s.out\n" % (job_directory, g, str(i*8), str((i+1)*8)))
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
                outdir = "%s/APS137_Archive/%s_OUT" % (os.getcwd(), 'ST-'+c[0]+'_ST-'+c[1])
                fh.writelines("cd %s\n" %outdir)
                fh.writelines("\n")
                fh.writelines("mcorr-xmfa-2clades %s/%s_OUT/MSA_%s_%s %s/%s_OUT/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                              % (wrkd, 'NGS_'+c[0], g, 'ST-'+c[0], wrkd, 'NGS_'+c[1], g, 'ST-'+c[1], outdir, 'ST-'+c[0]+'_ST-'+c[1], g))
                fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true\n"
                              % (outdir, 'ST-'+c[0]+'_ST-'+c[1], g, outdir, 'ST-'+c[0]+'_ST-'+c[1], g))
                fh.writelines("\n")
        os.system("sbatch %s" %job_file)
        print('submitted %s job #%s' %(g, i))
        time.sleep(2)
print('%s total submissions' %count)