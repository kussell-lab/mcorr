#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=SPpangenome
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=SPpangenomecollect_slurm%j.out

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
corrdir=${projectdir}/${projectdir}/APS160_SP_Archive/corethreshold95
outdir=${projectdir}/APS225_SP_Archive
cd ${outdir}

echo "let's rock"
collectFitAIC APS225_SP_pg_mcorr_results --corr_dir=${corrdir} --out_dir=${outdir}

