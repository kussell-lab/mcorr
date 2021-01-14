#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem=32GB
#SBATCH --job-name=APS160cluster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS160_cluster_slurm%j.out

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
jobdir=/scratch/aps376/recombo/APS160clusterseqs
strains=/scratch/aps376/recombo/APS160_SP_Archive/APS160_strain_list
outdir=/scratch/aps376/recombo/APS160_SP_Archive
##outfile
out="APS160_0113_SP_10thpercentile_cutoff"
##cutoff
cutoff=10
##minimum number of strains per cluster
min=100

cd ${outdir}

echo "let's rock"
clusterSequences APS160distancematrix.npy ${strains} ${out} --percentile=${cutoff} --min_size=${min}

