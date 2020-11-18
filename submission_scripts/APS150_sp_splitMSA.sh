#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=SP_splits
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_SP_splits.out


JOBDIR=/scratch/aps376/recombo/APS150_splitMSA
MSA=/scratch/aps376/recombo/APS150_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
ARCHIVE=/scratch/aps376/recombo/APS150_SP_Archive/
OUTDIR=${ARCHIVE}/SP_split_MSA
projectdir=/scratch/aps376/recombo

cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

mkdir -p $ARCHIVE
mkdir -p $OUTDIR

cd ${JOBDIR}

echo "let's rock"
python3 splitMSA.py $MSA MSA_SP_MASTER $OUTDIR 24

