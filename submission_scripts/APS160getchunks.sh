#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=50GB
#SBATCH --job-name=APS160chunks
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS160chunks_slurm%j.out

projectdir=/scratch/aps376/recombo
cd ${projectdir}

###activate venv and get other dependencies

module purge;
source venv/bin/activate;
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;
##load Go
module load go/1.15.2
#put gobin on path
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

## job directory
jobdir=/scratch/aps376/recombo/APS160mcorr-dm-chunks
strains=/scratch/aps376/recombo/APS160splitGenome/APS160_strain_list
##mcp output folder
mcp=/scratch/aps376/recombo/APS158_SP_Archive/0104_SP_mps_dists/
##cutoff
cutoff=10
##path to MSA
MSA=/scratch/aps376/recombo/APS158_SP_Archive/SP_MASTER_OUT/MSA_SP_PANGENOME_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS160_SP_Archive
cd ${outdir}

echo "let's rock"
mcorr-dm-chunks ${mcp} ${strains} APS160distancematrix --num-workers=16

