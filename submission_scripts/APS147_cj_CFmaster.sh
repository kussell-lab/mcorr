#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00
#SBATCH --mem=16GB
#SBATCH --job-name=APS147cjMSA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS147_CJ_CF_MSA_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

## job directory
jobdir=/scratch/aps376/recombo/APS147_avgrates

cd ${jobdir}

echo "let's rock"
python3 APS147_cj_CF_MSA.py

