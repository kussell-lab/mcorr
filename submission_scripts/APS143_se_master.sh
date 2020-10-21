#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=APS143se_MSA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS143se_MSA_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

## job directory
jobdir=/scratch/aps376/recombo/APS143_se_MSA

cd ${jobdir}

echo "let's rock"
python3 aps143_sentericaMSA.py