#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=APS148paMSA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS148_MSAsplits_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

## job directory
jobdir=/scratch/aps376/recombo/APS148_pa_split_MSA

cd ${jobdir}

echo "let's rock"
python3 APS148paSplitMSA.py

