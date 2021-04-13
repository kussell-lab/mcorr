#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=NMchunk
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=500chunks_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.15.7
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=$SCRATCH/recombo/APS183clusterseqs
##cutoff
splits=500
##outdir
outdir=${projectdir}/APS183_NM_Archive
list=${outdir}/strain_list
MSA=${outdir}/MSA_NM_MASTER_GAPFILTERED
cd ${outdir}

echo "let's rock"
chunkMSA ${MSA} ${list} ${splits} --chunk-folder="500chunks"

