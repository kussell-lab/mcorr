#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=3:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=HP_dm
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS192collectdm_slurm%j.out

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
jobdir=/scratch/aps376/recombo/APS192clusterseqs
archive=/scratch/aps376/recombo/APS192_HP_Archive
##mcp output
mcp=${archive}/APS192_HP_mcp_out.csv
strains=${archive}/strain_list
cd ${archive}
echo "let's rock"
mcorr-dm ${mcp} APS192distancematrix

