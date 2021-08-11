#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=APS226cluster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS226seqcluster_pt1_slurm%j.out

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
##archive
archive=/scratch/aps376/recombo/APS169_SP_Archive
##outdir
outdir=/scratch/aps376/recombo/APS226_SP_Archive
#strain list
strains=${archive}/strain_list
##outfile
out="APS226_SP_0.1thpercentile_cutoff"
##cutoff
cutoff=0.1

mkdir -p ${outdir}
cd ${outdir}
echo "let's rock"
clusterSequences ${archive}/APS169distancematrix.npy ${strains} ${out} --percentile=${cutoff} --min_size=100

