#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=APS156clusters
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS156clusters_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.13.6
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=/scratch/aps376/APS156makeclusters
##mcp output
mcp=/scratch/aps376/recombo/APS156_SP_Archive/1201_SP_mps_dists/APS156_201202_SP_all_dists.csv
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS156_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS156_SP_Archive
cd ${jobdir}

echo "let's rock"
makeClusters ${mcp} ${outdir} ${MSA}

