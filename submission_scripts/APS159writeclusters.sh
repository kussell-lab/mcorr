#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=APS159writeclusters
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS159writeclusters_slurm%j.out

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

##path to MSA
MSA=/scratch/aps376/recombo/APS159_NG_Archive/NG_MASTER_OUT/MSA_NG_PANGENOME_MASTER
##outdir
outdir=/scratch/aps376/recombo/APS159_NG_Archive
##cluster-dict
clusterlist=/scratch/aps376/recombo/APS159_NG_Archive/cluster_list

cd ${outdir}

echo "let's rock"
write-cluster-msa ${MSA} ${clusterlist}

