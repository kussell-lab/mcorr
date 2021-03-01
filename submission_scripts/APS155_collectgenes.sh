#!/bin/bash
#SBATCH --job-name=APS154collector
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --time=2:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu


##INPUTS
WRKD=${SCRATCH}/recombo/salmonella
#SRC is the directory your SRA accession files and genome files are within.
SRC=/home/aps376/salmonella
Archive=${SCRATCH}/recombo/APS155_SE_Archive/SE_MASTER_OUT
FASTA=${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.fna
GFF=${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.gff
sra_list=APS155_finalpiles

mkdir -p $Archive

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


##MSA stands for multi sequence alignment in the below
  #the '$1' command tells it to grab the argument of pipe_dream

echo "let's rock"

CollectGeneAlignments ${sra_list} ${GFF} ${WRKD} ${Archive}/MSA_SE_MASTER --appendix ".pileup.fasta" --progress

