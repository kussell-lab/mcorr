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
job_directory = "/scratch/aps376/APS137_MLST" #%os.getcwd()
ngs = "/scratch/aps376/ngonorrhoeae" #os.getcwd()
#scratch = os.environ['SCRATCH']
# mkdir_p(job_directory)
# mkdir_p(ngs)

'make all possible combos from sero_list file, and sort them alphabetically'
sra_list_file = '/scratch/aps376/NGS_SRR_Acc_List.txt'
#sra_list_file = '/Users/asherpreskasteinberg/go/src/github.com/apsteinberg/mcorr/core_flex/NGS_SRR_Acc_List.txt'

#read the list of serotypes
sra_list = []
with open(sra_list_file, 'r') as reader:
    for line in reader:
        sra_list.append(line.rstrip())

#for i in np.arange(0, 193):
for i in np.arange(0, 193):
    job_file = os.path.join(job_directory, "NGS_SRA_job%s.sh" % i)
    # outdir = ngs
    # mkdir_p(outdir)
    # os.system('cd %s' %outdir)
    accession_list = sra_list[(i*10):((i+1)*10)]

    with open(job_file, "w+") as fh:
        fh.writelines("#!/bin/bash\n")
        fh.writelines("#SBATCH --job-name=NGS_SRA_job%s\n" % i)
        fh.writelines("#SBATCH --nodes=1\n")
        fh.writelines("#SBATCH --cpus-per-task=8\n")
        fh.writelines("#SBATCH --time=2:00:00\n")
        fh.writelines("#SBATCH --mem=16GB\n")
        fh.writelines("#SBATCH --mail-type=END,FAIL\n")
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("#SBATCH --output=%s/NGS_MLST_job%s.out\n" % (job_directory, i))
        fh.writelines("#SBATCH --mail-user=aps376@nyu.edu\n")
        fh.writelines("\n")
        #load modules
        fh.writelines("module purge\n")
        fh.writelines("cd %s\n" %job_directory)
        fh.writelines("source venv/bin/activate;\n")
        fh.writelines("export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;\n")
        fh.writelines("\n")
        fh.writelines("cd %s\n" %ngs)
        fh.writelines("\n")
        for accession in accession_list:
            fh.writelines('stringMLST.py --predict -P neisseria/nmb -o '
                          + str(accession) +'_MLST -1 '+ str(accession) +'_1.fastq -2 '
                          + str(accession) +'_2.fastq\n')
    os.system("sbatch %s" %job_file)
    print('submitted NGS_MLST_%s' %i)
    time.sleep(2)