#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:30:00
#SBATCH --mem=8GB
#SBATCH --job-name=PA_splits
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_PA_splits.out


JOBDIR=/scratch/aps376/recombo/APS151_splitMSA
MSA=/scratch/aps376/recombo/APS148_PA_Archive/PA_MASTER_OUT/MSA_PA_MASTER
OUTDIR=/scratch/aps376/recombo/APS148_PA_Archive/PA_split_MSA
projectdir=/scratch/aps376/recombo

cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;


mkdir -p $OUTDIR

cd ${JOBDIR}

echo "let's rock"
python3 splitMSA.py $MSA MSA_PA_MASTER $OUTDIR 10

