#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --time=24:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=SRAs
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=SRA_slurm_%j.out

module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.7.3 ##do 3.7.3!
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
module load sra-tools/2.10.5
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

mkdir -p $SCRATCH/recombo/saureus

#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

python3 $HOME/ReferenceAlignmentGenerator/FetchSRA_APS2.py APS154_SRA_list $SCRATCH/recombo/saureus
