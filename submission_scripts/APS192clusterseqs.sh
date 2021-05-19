#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=HPcluster
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS192seqcluster_slurm%j.out

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
outdir=/scratch/aps376/recombo/APS192_HP_Archive
#strain list
strains=${outdir}/strains
##outfile
out="APS192_HP_10thpercentile_cutoff"
##cutoff
cutoff=10
min_strains=10

cd ${outdir}
echo "let's rock"
clusterSequences APS192distancematrix.npy ${strains} ${out} --percentile=${cutoff} --min_size=${min_strains}

