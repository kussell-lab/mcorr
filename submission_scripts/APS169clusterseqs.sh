#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=12:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=APS169cluster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS169seqcluster_slurm%j.out

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
##outdir
outdir=/scratch/aps376/recombo/APS169_SP_Archive
#strain list
strains=${outdir}/strain_list
##outfile
out="APS169_SP_10thpercentile_cutoff"
##cutoff
cutoff=10

cd ${outdir}
echo "let's rock"
clusterSequences APS169distancematrix.npy ${strains} ${out} --percentile=${cutoff} --min_size=100

