#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=APS150clusters
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS148_MSAsplits_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

## job directory
jobdir=/scratch/aps376/recombo/APS150_sp_clusters
##mcp output
mcp=${jobdir}/APS150_201106_SP_all_dists.csv
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS150_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS150_SP_Archive/cluster_MSAs
cd ${jobdir}

echo "let's rock"
python3 makeSeqClusters.py ${mcp} ${cutoff} ${MSA} ${outdir}

