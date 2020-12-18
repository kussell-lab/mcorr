#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=1:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=APS156geneMSA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS156mcorr_slurm%j.out

module purge
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

#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

## job directory
jobdir=/scratch/aps376/recombo/APS156geneMSA
strains=${jobdir}/APS156_finalpile
##mcp output
mcp=/scratch/aps376/recombo/APS156_SP_Archive/1201_SP_mps_dists/
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS156_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS156_SP_Archive/genes

echo "let's rock"
geneMSA ${MSA} --outdir=${outdir}
