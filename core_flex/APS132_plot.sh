#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=APS132_plot
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm_%j.out ##can be changed to "$OUTPUT" in the iterated code

module purge

module load python/intel/2.7.12 ##do 3.7.3!
module load pandas/intel/py2.7/0.20.3
module load seaborn/0.7.1

seros='sero_list_7'
filedir='/scratch/aps376/Archive/'

python APS132_get_divergence.py ${seros} ${filedir}
