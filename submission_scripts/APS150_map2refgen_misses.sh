#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=map2refgen
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_map2refgen.out

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

WRKDIR=$SCRATCH/recombo/spneumoniae
OUTDIR=$SCRATCH/recombo/APS150_SP_Archive
SRC=/home/aps376/APS150_spneumoniae
FASTA=$SRC/Reference/GCF_000007045.1_ASM704v1_genomic.fna
GFF=$SRC/Reference/GCF_000007045.1_ASM704v1_genomic.gff
LISTS=/scratch/aps376/recombo/APS150_map_misses/APS150_lastpiles

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

cd $WRKDIR

echo "let's rock"
bash MapRead2Reference.sh ${LISTS} $WRKDIR $FASTA


