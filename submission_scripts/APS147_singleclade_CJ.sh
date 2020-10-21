#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=APS147CJ
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS147CJ_slurm%j.out

module purge
module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
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

#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

OUTDIR=$SCRATCH/recombo/APS147_avgrates

echo "let's rock"
mcorr-xmfa ${OUTDIR}/MSA_CJ_CORE_Master_Sorted ${OUTDIR}/CJ_CORE_XMFA_OUT --show-progress &&
mcorr-fit ${OUTDIR}/CJ_CORE_XMFA_OUT.csv ${OUTDIR}/CJ_CORE_FIT_OUT || true

mcorr-xmfa ${OUTDIR}/MSA_CJ_FLEX_Master_Sorted ${OUTDIR}/CJ_FLEX_XMFA_OUT --show-progress &&
mcorr-fit ${OUTDIR}/CJ_FLEX_XMFA_OUT.csv ${OUTDIR}/CJ_FLEX_FIT_OUT || true

