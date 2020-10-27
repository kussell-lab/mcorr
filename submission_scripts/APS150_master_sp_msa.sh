#!/bin/bash
#SBATCH --job-name=APS148MSA
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=72:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu


##INPUTS
#The directory you want create serotype folders in. You need lots of space.
WRKD=${BEEGFS}
#SRC is the directory your SRA accession files and genome files are within.
SRC=/home/aps376/APS150_spneumoniae
Archive=${SCRATCH}/recombo/APS150_SP_Archive
FASTA=$SRC/Reference/GCF_000007045.1_ASM704v1_genomic.fna
GFF=$SRC/Reference/GCF_000007045.1_ASM704v1_genomic.gff
sra_list=${SRC}/APS150_downloaded_sra_list

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
mkdir ${Archive}/SP_MASTER_OUT
cd ${Archive}/SP_MASTER_OUT

echo "let's rock"

ReferenceAlignmentGenerate ${sra_list} ${WRKD}/spneumoniae ${FASTA} ${GFF} ${Archive}/SP_MASTER_OUT/MSA_SP_MASTER


