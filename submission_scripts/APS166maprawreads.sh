#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=0:30:00
#SBATCH --mem=4GB
#SBATCH --job-name=APS166bowtie
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=APS166bowtie_slurm%j.out

##INPUTS
projdir=/scratch/aps376/recombo
jobdir=/scratch/aps376/recombo/APS166rawreads
WRKD=/scratch/aps376/recombo/APS166_BS_Archive
prokka=${WRKD}/prokka
list=${WRKD}/assembly_accession_list
genomes=${WRKD}/genomes
roary_dir=${WRKD}/roary
ncpu=30

echo "Loading modules."
module load go/1.15.2
module load python/intel/3.8.6
module load parallel/20201022
module load samtools/intel/1.11
module load singularity/3.6.4
##load virtual environment
#cd ${projdir}
#source venv/bin/activate;
### maybe this line is the issue?
#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK;

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path
echo "Making everything in path."
#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

##set perl language variable; this will give you fewer warnings
export LC_ALL=C

##MSA stands for multi sequence alignment in the below
#the '$1' command tells it to grab the argument of pipe_dream

echo "let's rock"
cd ${WRKD}
logfile="roary.log"
roary='singularity exec /scratch/work/public/singularity/roary-3.13.0.sif roary'
#mkdir -p ${roary_dir}
${roary} -v -p ${ncpu} -f ${roary_dir} ${roary_dir}_gffs/*.gff &>$logfile
mv $logfile ${roary_dir}
echo "Completed running Roary, see log file in ${roary_dir}/roary.log"