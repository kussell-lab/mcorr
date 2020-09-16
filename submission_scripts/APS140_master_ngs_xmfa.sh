#!/bin/bash
#SBATCH --job-name=APS140MSA
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
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
SRC=/home/aps376/APS137_ngs
Archive='/scratch/aps376/APS140_Archive'
FASTA=$SRC/Reference/GCF_000020105.1_ASM2010v1_genomic.fna
GFF=$SRC/Reference/GCF_000020105.1_ASM2010v1_genomic.gff
sra_list='/scratch/aps376/NGS_SRR_Acc_List.txt'

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

mkdir -p ${Archive}

##MSA stands for multi sequence alignment in the below
  #the '$1' command tells it to grab the argument of pipe_dream
mkdir ${Archive}/NGS_MASTER_OUT
cd ${Archive}/NGS_MASTER_OUT

ReferenceAlignmentGenerate ${sra_list} ${WRKD}/ngonorrhoeae ${FASTA} ${GFF} ${Archive}/NGS_MASTER_OUT/MSA_NGS_MASTER

echo "let's rock"
