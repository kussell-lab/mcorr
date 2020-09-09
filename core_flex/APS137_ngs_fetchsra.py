#!/usr/bin/env python
import os
import time
from itertools import combinations
import numpy as np
from tqdm import tqdm


def mkdir_p(dir):
    'make a directory if doesnt exist'
    if not os.path.exists(dir):
        os.mkdir(dir)

#define directories
job_directory = "%s/APS137_fetchsra" %os.getcwd()
ngs = "%s/ngonorrhoeae" %os.getcwd()
#scratch = os.environ['SCRATCH']
mkdir_p(job_directory)
mkdir_p(ngs)

'make all possible combos from sero_list file, and sort them alphabetically'
sra_list_file = '/scratch/aps376/NGS_SRR_Acc_List.txt'
#sra_list_file = '/Users/asherpreskasteinberg/go/src/github.com/apsteinberg/mcorr/core_flex/NGS_SRR_Acc_List.txt'

#read the list of serotypes
sra_list = []
with open(sra_list_file, 'r') as reader:
    for line in reader:
        sra_list.append(line.rstrip())

for i in np.arange(0, 193):
    job_file = os.path.join(job_directory, "NGS_SRA_job%s.sh" % i)
    outdir = ngs
    mkdir_p(outdir)
    os.system('cd %s' %outdir)
    accession_list = sra_list[(i*10):((i+1)*10)]

    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=NGS_SRA_job%s\n" % i)
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=8\n")
        fh.writelines("#SBATCH --time=1:00:00\n")
        fh.writelines("#SBATCH --mem=4GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/NGS_SRA_job%s.out\n" % (job_directory, i))
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
        for accession in accession_list:
            fh.writelines('fasterq-dump ' + str(accession) + ' -O ' + outdir + ' -t ' + '$SCRATCH\n')
    os.system("sbatch %s" %job_file)
    print('submitted NGS_SRA_job_%s' %i)
    time.sleep(2)