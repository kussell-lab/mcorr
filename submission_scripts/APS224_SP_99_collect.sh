#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=SP99collect
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=SP99collect_slurm%j.out

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


##corrdir
corrdir=${projectdir}/APS169_SP_Archive/corethreshold99
outdir=${projectdir}/APS224_SP_99_Archive
cd ${outdir}

echo "let's rock"
collectFitAIC APS224_SP99_mcorr_results --corr_dir=${corrdir} --out_dir=${outdir}

