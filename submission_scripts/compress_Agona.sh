#!/bin/bash
#SBATCH --job-name=APS135MSA
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu

####
#written by Asher, 200731
#####

#The directory you want create serotype folders in. You need lots of space.
WRKD=/scratch/aps376

tar zcvf Agona.tar.gz Agona

