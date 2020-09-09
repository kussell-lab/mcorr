#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=getST
#SBATCH --cpus-per-task=4
#SBATCH --mem=4GB
#SBATCH --time=0:30:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm_%j.out ##can be changed to "$OUTPUT" in the iterated code

jobdir='/scratch/aps376/APS137_MLST'

cd jobdir
module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

python3 ${jobdir}/APS137_getST.py