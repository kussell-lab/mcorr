#!/bin/bash
#SBATCH --job-name=APS135MSA
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --array=1-18
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

###
#BACKGROUND
#This script makes gene aligment files
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
SRC=/home/aps376/APS135_salmonella
Archive='/scratch/aps376/APS135_Archive'

echo "Loading modules."
module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.6.3 ##do 3.7.3!
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

mkdir ${Archive}

##MSA stands for multi sequence alignment in the below
pipe_dream () {
  #the '$1' command tells it to grab the argument of pipe_dream
	SERO=$1
	echo files_$SERO
	mkdir ${Archive}/${SERO}_OUT
	cd ${Archive}/${SERO}_OUT  
	IFS='-'
	read -a strarr <<< "$1"
	echo "main sero: ${strarr[0]}"
	#cd "/scratch/aps376/Archive/${strarr[0]}_out"

        CollectGeneAlignments ${SRC}/SRA_files/sra_accession_$sero ${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.gff ${WRKD}/${strarr[0]} ${Archive}/${SERO}_OUT/MSA_${SERO} --appendix ".pileup.fasta" --progress
}

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#this line grabs the line from the list of serotypes in sero_list
#for example when SLURM_ARRAY_TASK_ID=1, you grab the first line
sero=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sero_list)
echo "Running Serotype: " $sero
pipe_dream $sero
