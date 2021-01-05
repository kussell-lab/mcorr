#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=APS158cluster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS158cluster_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.15.2
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=/scratch/aps376/recombo/APS158makeclusters
strains=/scratch/aps376/recombo/APS158_spneumoniae/APS156_completepiles
##mcp output folder
mcp=/scratch/aps376/recombo/APS158_SP_Archive/0104_SP_mps_dists/
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS158_SP_Archive/SP_MASTER_OUT/MSA_SP_PANGENOME_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS158_SP_Archive
cd ${outdir}

echo "let's rock"
clusterSequences APS158distancematrix.npy ${strains}

