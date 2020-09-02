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
job_directory = "%s/APS135_v1" %os.getcwd()
archive = "%s/APS135_Archive" %os.getcwd()
#scratch = os.environ['SCRATCH']
mkdir_p(job_directory)
mkdir_p(archive)

# combolist = [('Kentucky','Oranienburg','FLEX'),
#              ('Newport', 'Virchow', 'FLEX'),
#              ('Agona','Typhimurium','FLEX'),
#              ('Braenderup','Typhimurium','FLEX'),
#              ('Infantis','Typhimurium','FLEX'),
#              ('Agona','Enteritidis','FLEX'),
#              ('Kentucky','Typhimurium','FLEX'),
#              ('Braenderup','Enteritidis','FLEX')]

# combolist = [('Enteritidis','Typhimurium','FLEX'),
#              ('Newport','Typhimurium','FLEX'),
#              ('Oranienburg','Typhimurium','FLEX'),
#              ('Agona','Oranienburg','FLEX'),
#              ('Braenderup','Oranienburg','FLEX'),
#              ('Enteritidis','Infantis','FLEX'),
#              ('Enteritidis','Kentucky','FLEX'),
#              ('Enteritidis','Newport','FLEX'),
#              ('Enteritidis','Stanley','FLEX'),
#              ('Kentucky','Oranienburg','FLEX')
#              ]
combolist = [('Agona-EBG54','Newport-EBG3','FLEX'),
             ('Agona-EBG54','Oranienburg-EBG41','FLEX'),
             ('Agona-EBG54','Stanley-EBG29','FLEX'),
             ('Agona-EBG54','Virchow-EBG70','FLEX'),
             ('Braenderup-EBG24','Kentucky-EBG56','FLEX'),
             ('Agona-EBG54','Kentucky-EBG56','FLEX'),
             ('Agona-EBG54','Virchow-EBG9','FLEX'),
             ('Braenderup-EBG24','Virchow-EBG70','FLEX'),
             ('Infantis-EBG31','Kentucky-EBG164','FLEX'),
             ('Agona-EBG54','Oranienburg-EBG44','FLEX'),
             ('Infantis-EBG31','Oranienburg-EBG41','FLEX'),
             ('Braenderup-EBG24','Stanley-EBG29','FLEX'),
             ('Braenderup-EBG24','Typhi-EBG13','FLEX'),
             ('Kentucky-EBG164','Kentucky-EBG56','FLEX'),
             ('Braenderup-EBG24','Newport-EBG7','FLEX'),
             ('Kentucky-EBG164','Oranienburg-EBG44','FLEX'),
             ('Kentucky-EBG164','Oranienburg-EBG41','FLEX'),
             ('Kentucky-EBG56','Virchow-EBG70','FLEX'),
             ('Newport-EBG2','Newport-EBG3','FLEX'),
             ('Newport-EBG2','Oranienburg-EBG41','FLEX'),
             ('Enteritidis-EBG4','Virchow-EBG9','FLEX'),
             ('Newport-EBG3','Oranienburg-EBG44','FLEX'),
             ('Newport-EBG3','Oranienburg-EBG41','FLEX'),
             ('Kentucky-EBG164','Stanley-EBG29','FLEX'),
             ('Infantis-EBG31','Newport-EBG7','FLEX'),
             ('Kentucky-EBG164','Typhimurium-EBG1','FLEX'),
             ('Kentucky-EBG164','Virchow-EBG9','FLEX'),
             ('Enteritidis-EBG4','Stanley-EBG29','FLEX'),
             ('Oranienburg-EBG41','Virchow-EBG9','FLEX'),
             ('Enteritidis-EBG4','Infantis-EBG31','FLEX'),
             ('Oranienburg-EBG44','Virchow-EBG9','FLEX'),
             ('Enteritidis-EBG4','Newport-EBG3','FLEX'),
             ('Newport-EBG2','Typhimurium-EBG138','FLEX'),
             ('Newport-EBG7','Stanley-EBG29','FLEX'),
             ('Typhi-EBG13','Typhimurium-EBG138','FLEX'),
             ('Agona-EBG54','Enteritidis-EBG4','CORE'),
             ('Newport-EBG7','Typhimurium-EBG1','FLEX'),
             ('Braenderup-EBG24','Enteritidis-EBG4','CORE'),
             ('Braenderup-EBG24','Typhimurium-EBG1','CORE'),
             ('Enteritidis-EBG4','Typhimurium-EBG1','CORE'),
             ('Kentucky-EBG56','Typhimurium-EBG1','CORE')]
wrkd = '/scratch/aps376/APS135_Archive'

gene = ['CORE', 'FLEX']

for g in gene:
    for c in combolist:
        if c[2] != g:
            continue
        job_file = os.path.join(job_directory, "%s_%s.sh" % (g, c[0]+c[1]))
        outdir = "%s/APS135_Archive/%s_OUT" % (os.getcwd(), c[0]+c[1])
        mkdir_p(outdir)
        os.system('cd %s' %outdir)

        with open(job_file, "w+") as fh:
            fh.writelines("#!/bin/bash\n")
            fh.writelines("#SBATCH --job-name=%s_%s\n" % (g, c[0]+c[1]))
            fh.writelines("#SBATCH --nodes=1\n")
            fh.writelines("#SBATCH --cpus-per-task=16\n")
            fh.writelines("#SBATCH --time=24:00:00\n")
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
            fh.writelines("mcorr-xmfa-2clades %s/%s_OUT/MSA_%s_%s %s/%s_OUT/MSA_%s_%s %s/%s_%s_XMFA_OUT &&\n"
                          % (wrkd, c[0], g, c[0], wrkd, c[1], g, c[1], outdir, c[0]+c[1], g))
            fh.writelines("mcorr-fit %s/%s_%s_XMFA_OUT.csv %s/%s_%s_FIT_OUT || true"
                          % (outdir, c[0]+c[1], g, outdir, c[0]+c[1], g))
        os.system("sbatch %s" %job_file)
        print('submitted %s for %s' %(g, c[0]+c[1]))
        time.sleep(2)