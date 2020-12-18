#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00
#SBATCH --mem=8GB
#SBATCH --job-name=APS154chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS150clusters_slurm%j.out

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
jobdir=$BEEGFS/APS154_chunkMSA
##cutoff
splits=200
##path to MSA
MSA=/beegfs/aps376/APS154_SA_Archive/SA_MASTER_OUT/MSA_SA_MASTER
##outdir
outdir=/beegfs/aps376/APS154_SA_Archive/
cd ${outdir}

echo "let's rock"
chunkMSA ${MSA} ${splits}

