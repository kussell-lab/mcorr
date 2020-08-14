#!/bin/bash
#SBATCH --job-name=FIT_CORE_FLEX
#SBATCH --cpus-per-task=8
#SBATCH --mem=32GB
#SBATCH --array=1-8
#SBATCH --time=6:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

###
#BACKGROUND
#This script assumes you've already run the normal 'pipeline' script and the 'RefGen_Core_Flex_Split.py' script.

#NOTES
#Minor untested edits have been made to this file to make it more readable.
#If something breaks, then contact me for the original or assistance. -Spencer
##with edits by Asher
###
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

##INPUTS
#The directory you want create serotype folders in. You need lots of space.
WRKD=/scratch/aps376

pipe_dream () {
	SERO=$1
	cd ${WRKD}/Archive/${SERO}_OUT

	for gt in {'CORE','FLEX'}
	do
	echo $gt
        mcorr-xmfa ${WRKD}/Archive/${SERO}_OUT/REFGEN_${gt}_${SERO} ${WRKD}/Archive/${SERO}_OUT/${SERO}_${gt}_XMFA_OUT &&
        mcorr-fit ${WRKD}/Archive/${SERO}_OUT/${SERO}_${gt}_XMFA_OUT.csv ${WRKD}/Archive/${SERO}_OUT/${SERO}_${gt}_FIT_OUT || true #Needs X11 forwarding. Can/should fix this.
	done
}

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
sero=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sero_list)
echo "Running Serotype: " $sero
pipe_dream $sero
