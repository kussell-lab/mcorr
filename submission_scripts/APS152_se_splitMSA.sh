#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=SE_splits
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_SE_splits.out


JOBDIR=/scratch/aps376/recombo/APS152_splitMSA
MSA=/scratch/aps376/recombo/APS143_1008_senterica_Archive/MSA_Master_Sorted
ARCHIVE=/scratch/aps376/recombo/APS152_se_Archive/
OUTDIR=/scratch/aps376/recombo/APS152_se_Archive/se_split_MSA
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
python3 splitMSA.py $MSA MSA_Master_Sorted $OUTDIR 250

