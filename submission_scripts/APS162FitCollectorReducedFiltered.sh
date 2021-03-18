#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:30:00
#SBATCH --mem=4GB
#SBATCH --job-name=APS162fitcollector
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS162fitcollector_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.15.7
module load python/intel/3.8.6
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=/scratch/aps376/recombo/APS160.5fitcollector

##outdir
archive=/scratch/aps376/recombo/APS162_SP_Archive
##outdir
outdir=${archive}/APS156filteredreduced_threshold90

cd ${outdir}

echo "let's rock"
FitCollector ${outdir} --lmfitSuffix="_FIT_OUT_lmfit_report.csv"

