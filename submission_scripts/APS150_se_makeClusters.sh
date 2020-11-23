#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=APS150clusters
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS150clusters_slurm%j.out

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
jobdir=/beegfs/aps376/APS152_se_go_clusters
##mcp output
mcp=/beegfs/aps376/APS152_se_Archive/APS152_201120_SE_all_dists.csv
##cutoff
cutoff=10
##path to MSA
MSA=/beegfs/aps376/APS143_1008_senterica_Archive/MSA_Master_Sorted
##outdir
outdir=/beegfs/aps376/APS152_se_Archive
cd ${jobdir}

echo "let's rock"
makeClusters ${mcp} ${outdir} ${MSA}

