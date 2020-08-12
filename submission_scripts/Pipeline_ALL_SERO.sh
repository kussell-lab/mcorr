#!/bin/bash
#SBATCH --job-name=Corr_Profiles
#SBATCH --cpus-per-task=8
#SBATCH --mem=32000
#SBATCH --array=1-10
#SBATCH -t2-0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

###
#BACKGROUND
#This script runs the entire mcorr pipeline using 'ReferenceAlignmentGenerate'
#This script assumes you have one folder containing two folders and one file:
#-SRA_files, a folder containing SRA accession files
#-Reference, a folder containing the gff and fna files for your genome
#-sero_list, a file where each line is the name of a serotype. The names should match the
# names within the SRA_files folder. The batch array length corresponds to the number of lines
# you would like to read from this file.
#NOTES
#Minor edits have been added to this file to make the file paths more readable.
#If something breaks, then be suspicious of the file paths.
###

##INPUTS
#The directory you want create serotype folders in. You need lots of space.
WRKD=/scratch/aps376
#SRC is the directory your SRA accession files and genome files are within.
SRC=/home/aps376/salmonella

echo "Loading modules."
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
echo "Making everything in path."
#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

mkdir Archive

pipe_dream () {
  #the '$1' command tells it to grab the argument of pipe_dream
	SERO=$1
	echo files_$SERO
	mkdir ${WRKD}/$SERO
	mkdir ${WRKD}/Archive/${SERO}_OUT
	cd ${WRKD}/Archive/${SERO}_OUT

        ReferenceAlignmentGenerate ${SRC}/SRA_files/sra_accession_$sero ${WRKD}/$sero ${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.fna ${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.gff ${WRKD}/Archive/${SERO}_OUT/REFGEN_$SERO &&
        mcorr-xmfa ${WRKD}/Archive/${SERO}_OUT/REFGEN_$SERO ${WRKD}/Archive/${SERO}_OUT/${SERO}_XMFA_OUT &&
        mcorr-fit ${WRKD}/Archive/${SERO}_OUT/${SERO}_XMFA_OUT.csv ${WRKD}/Archive/${SERO}_OUT/${SERO}_FIT_OUT || true #Needs X11 forwarding. Can/should fix this.
}

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#this line grabs the line from the list of serotypes in sero_list
#for example when SLURM_ARRAY_TASK_ID=1, you grab the first line
sero=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sero_list)
echo "Running Serotype: " $sero
pipe_dream $sero
