#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=20GB
#SBATCH --job-name=APS143se_MSA
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS143se_MSA_slurm%j.out

module load python3/intel/3.6.3 ##do 3.7.3!

echo "let's rock"
python3 aps143_sentericaMSA.py