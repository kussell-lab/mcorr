#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=0:30:00
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
jobdir=/scratch/aps376/recombo/APS153_sa_go_clusters
##mcp output
mcp=/scratch/aps376/recombo/APS141_1008_Archive/saureus_mp_MASTER_XMFA_OUT.csv
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS141_1008_Archive/APS141_saureus_MASTER_MSA
##outdir
outdir=/scratch/aps376/recombo/APS153_saureus_Archive
cd ${jobdir}

echo "let's rock"
makeClusters ${mcp} ${outdir} ${MSA}

