#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:30:00
#SBATCH --mem=32GB
#SBATCH --job-name=APS143pairs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS143pairs_slurm%j.out


projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

## job directory
jobdir=/scratch/aps376/recombo/APS143_se_pairs

cd ${jobdir}

echo "let's rock"
python3 getIsolatePairs.py