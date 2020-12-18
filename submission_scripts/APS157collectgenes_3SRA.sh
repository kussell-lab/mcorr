#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=APS157collectgenes
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS157genes_3SRA_slurm%j.out


##INPUTS
jobdir=/scratch/aps376/recombo/APS157map2multirefgen
fastq_list=${jobdir}/APS157_genomestats/201218_SRA_list
GFF_list=${jobdir}/APS157_genomestats/refgenome_gff_list
GFF_dir=/scratch/aps376/recombo/APS157_SP_Archive/genome_assemblies_genome_gff/ncbi-genomes-2020-12-15
NCBIREF=GCF_000007045.1_ASM704v1_genomic.gff
WRKDIR=/scratch/aps376/recombo/APS157_spneumoniae
OUTDIR=/scratch/aps376/recombo/APS157_SP_Archive
OUT=${OUTDIR}/APS157_3SRA_XMFA
#WRKRS=14

mkdir -p ${OUTDIR}

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

CollectMultiGeneAlignments ${fastq_list} ${GFF_list} ${GFF_dir} ${NCBIREF} ${WRKDIR} ${OUT} --progress --appendix=".pileup.fasta"

