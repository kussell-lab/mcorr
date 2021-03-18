#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --mem=100GB
#SBATCH --job-name=APS170chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS170collectchunks_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.15.7
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=/scratch/aps376/recombo/APS170clusterseqs
archive=/scratch/aps376/recombo/APS170_SA_Archive
##mcp output
mcp=${archive}/0309_SA_mps_dists/
strains=${archive}/strain_list
cd ${archive}
echo "let's rock"
mcorr-dm-chunks ${mcp} ${strains} APS170distancematrix

