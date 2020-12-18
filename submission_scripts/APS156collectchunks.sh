#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=SP_chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_SE_splits.out


JOBDIR=/scratch/aps376/recombo/APS156collectchunks
projectdir=/scratch/aps376/recombo

cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

cd ${JOBDIR}

echo "let's rock"
python3 APS156collectChunks.py

