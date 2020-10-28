#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=PA_MPS
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_PA_MPS.out

module load git/gnu/2.16.2
module load go/1.13.6 #try go/1.13.6
module load python3/intel/3.6.3 ##do 3.7.3!
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
#module load sra-tools/intel/2.9.6 #try 2.9.6
module load sra-tools/2.10.5
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

MSA1=/scratch/aps376/recombo/APS148_PA_Archive/PA_MASTER_OUT/MSA_PA_MASTER
MSA2=/scratch/aps376/recombo/APS148_PA_Archive/PA_split_MSA/MSA_PA_MASTER_split0
OUTDIR=/scratch/aps376/recombo/APS148_PA_Archive/PA_mps_dists

mkdir -p $OUTDIR

#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

cd $OUTDIR

echo "let's rock"
mcorr-pair-sync $MSA1 $MSA2 $OUTDIR/PA_MPS_0_MASTER_XMFA_OUT.csv  --max-corr-length=3

