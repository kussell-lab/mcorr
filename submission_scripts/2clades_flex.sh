#!/bin/bash
#SBATCH --job-name=2cladeflex
#SBATCH --cpus-per-task=4
#SBATCH --mem=64GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

####
#written by Asher, 200731
#####
echo "Loading modules."
module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.7.3 ##might need to roll back to 3.6.3? unclear
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
module load sra-tools/2.10.5
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path
echo "Making everything in path."
#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

#The directory you want create serotype folders in. You need lots of space.
WRKD=/scratch/aps376

mkdir ${WRKD}/Archive/KentuckyNewport_OUT
OUTDIR=${WRKD}/Archive/KentuckyNewport_OUT

cd ${OUTDIR}

mcorr-xmfa-2clades ${WRKD}/Archive/Kentucky_OUT/REFGEN_FLEX_Kentucky ${WRKD}/Archive/Newport_OUT/REFGEN_FLEX_Newport ${OUTDIR}/KentuckyNewport_FLEX_XMFA_OUT &&
echo "time to FLEX"
mcorr-fit ${OUTDIR}/KentuckyNewport_FLEX_XMFA_OUT.csv ${OUTDIR}/KentuckyNewport_FLEX_FIT_OUT || true

echo "The party's over, go home"
