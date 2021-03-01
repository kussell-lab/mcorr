#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00
#SBATCH --mem=200GB
#SBATCH --job-name=mcorr_SC2
#SBATCH --mail-type=END
##SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=mcorrall_slurm_%j.out

# Load modules
module purge
module load samtools/intel/1.11
module load sra-tools/2.10.9
module load parallel/20201022
module load python/intel/3.8.6
module load smalt/intel/0.7.6
module load go/1.15.7
module load bowtie2/2.4.2
module load bedtools/intel/2.29.2

projdir=/scratch/aps376/recombo
archive=${projdir}/APS168_SC2_Archive
MSA=${archive}/MSA_SC2_MASTER

cd ${projdir}
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

#put go on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

mcorr-xmfa ${MSA} ${archive}/SC2_ALLSEQS_XMFA_OUT --show-progress &&
#mcorr-fit ${archive}/SC2_ALLSEQS_XMFA_OUT.csv ${archive}/SC2_ALLSEQS_FIT_OUT || true