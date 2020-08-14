#!/bin/bash
#SBATCH --job-name=Corr_Profiles
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH -t2-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

#The directory you want create serotype folders in. You need lots of space.
WRKD=/scratch/aps376

module load sra-tools/2.10.5

for sra in 'SRR1968826' 'SRR1968929' 'SRR1965464' 'SRR1966755' 'SRR1967009' 'SRR3049525'
do
  fasterq-dump ${sra} -O $SCRATCH/Enteritidis -t $SCRATCH
done