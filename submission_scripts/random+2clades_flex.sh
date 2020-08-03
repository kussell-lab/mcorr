#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=2cladetests
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=2clades_slurm_%j.out ##can be changed to "$OUTPUT" in the iterated code

####
#written by Asher, 200731
#####
echo "Loading modules."
module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.6.3 ##might need to roll back to 3.6.3? unclear
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
SRC=/home/aps376/recombo
Newport=/scratch/aps376/Archive/Newport_OUT/
Kentucky=/scratch/aps376/Archive/Kentucky_OUT/

mkdir ${WRKD}/Archive/betweencladestest
OUTDIR=${WRKD}/Archive/betweencladestest

cd ${Kentucky}
python3 ${SRC}/random_choose_strains.py REFGEN_Kentucky Kentucky_50.xmfa 50
python3 ${SRC}/random_choose_strains.py REFGEN_CORE_Kentucky Kentucky_CORE_50.xmfa 50
python3 ${SRC}/random_choose_strains.py REFGEN_FLEX_Kentucky Kentucky_FLEX_50.xmfa 50

cd ${Newport}
python3 ${SRC}/random_choose_strains.py REFGEN_Newport Newport_50.xmfa 50
python3 ${SRC}/random_choose_strains.py REFGEN_CORE_Newport Newport_CORE_50.xmfa 50
python3 ${SRC}/random_choose_strains.py REFGEN_FLEX_Newport Newport_FLEX_50.xmfa 50

cd ${OUTDIR}
mcorr-xmfa ${Kentucky}/Kentucky_50.xmfa ${OUTDIR}/Kentucky_50_XMFA_OUT &&
mcorr-fit ${OUTDIR}/Kentucky_50_XMFA_OUT.csv ${OUTDIR}/Kentucky_50_FIT_OUT || true
echo "still going"
mcorr-xmfa ${Newport}/Newport_50.xmfa ${OUTDIR}/Newport_50_XMFA_OUT &&
mcorr-fit ${OUTDIR}/Newport_50_XMFA_OUT.csv ${OUTDIR}/Newport_50_FIT_OUT || true
echo "not dead yet"
mcorr-xmfa-2clades ${Kentucky}/Kentucky_FLEX_50.xmfa ${Newport}/Newport_FLEX_50.xmfa ${OUTDIR}/KN_FLEX_XMFA_OUT &&
mcorr-fit ${OUTDIR}/KN_FLEX_XMFA_OUT.csv ${OUTDIR}/KN_FLEX_FIT_OUT || true
echo "still going"
mcorr-xmfa-2clades ${Kentucky}/Kentucky_CORE_50.xmfa ${Newport}/Newport_CORE_50.xmfa ${OUTDIR}/KN_CORE_XMFA_OUT &&
mcorr-fit ${OUTDIR}/KN_CORE_XMFA_OUT.csv ${OUTDIR}/KN_FLEX_FIT_OUT || true

echo "The party's over, go home"
