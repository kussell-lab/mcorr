#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=250chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=chunks_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.13.6
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=$SCRATCH/recombo/APS156chunkMSA
##cutoff
splits=250
##path to MSA
MSA=/scratch/aps376/recombo/APS156_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS156_SP_Archive/
cd ${outdir}

echo "let's rock"
chunkMSA ${MSA} ${splits} --chunk-folder="250chunks"

