#!/bin/bash
#SBATCH --nodes=1
#SBATCH --job-name=CF_split
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=3:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm_%j.out ##can be changed to "$OUTPUT" in the iterated code

SRC=/home/aps376/recombo

module purge
module load python3/intel/3.6.3 ##do 3.7.3!

python3 $SRC/RefGen_Core_Flex_Split_APS1.py