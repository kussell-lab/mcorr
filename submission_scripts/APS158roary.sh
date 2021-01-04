#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=APS158roary
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS158roary_slurm%j.out

##INPUTS
jobdir=/scratch/aps376/recombo/APS158map2multirefgen
WRKDIR=/scratch/aps376/recombo/APS158_SP_Archive

echo "Loading modules."
module load git/gnu/2.16.2
module load go/1.15.2          #try go/1.13.6
module load python/intel/3.8.6 ##do 3.7.3!
module load parallel/20201022
module load prokka/1.12
module load muscle/intel/3.8.31
module load sra-tools/2.10.5
module load samtools/intel/1.11
module load smalt/intel/0.7.6
module load singularity/3.6.4

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

singularity exec ${HOME}/roary_latest.sif roary -p 8 -v -f ${WRKDIR}/1223_roary ${WRKDIR}/1223_roary_gffs/*.gff -e -n
