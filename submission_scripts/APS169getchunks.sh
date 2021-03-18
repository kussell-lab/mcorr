#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=1:00:00
#SBATCH --mem=100GB
#SBATCH --job-name=APS169chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS169clusters_slurm%j.out

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
jobdir=/scratch/aps376/recombo/APS169clusterseqs
##mcp output
mcp=/scratch/aps376/recombo/APS169_SP_Archive/0309_SP_mps_dists/
##outdir
outdir=/scratch/aps376/recombo/APS169_SP_Archive
strains=${outdir}/strain_list
cd ${outdir}
echo "let's rock"
mcorr-dm-chunks ${mcp} ${strains} APS169distancematrix

