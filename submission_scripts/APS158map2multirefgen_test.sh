#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=APS158map
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS158map_slurm%j.out


##INPUTS
jobdir=/scratch/aps376/recombo/APS158map2multirefgen
fastq_list=${jobdir}/APS157_genomestats/201221_SRA_test
FASTA_list=${jobdir}/APS157_genomestats/refgenome_fasta_list
FASTA_dir=/scratch/aps376/recombo/APS157_SP_Archive/APS157_SP_genome_assemblies_genome_fasta/ncbi-genomes-2020-12-15
WRKDIR=/scratch/aps376/recombo/APS158_sp_test
WRKRS=14

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

##set perl language variable; this will give you fewer warnings
export LC_ALL=C


##MSA stands for multi sequence alignment in the below
  #the '$1' command tells it to grab the argument of pipe_dream

echo "let's rock"

map2multirefgen ${fastq_list} ${FASTA_list} ${WRKDIR} --workers=${WRKRS}

