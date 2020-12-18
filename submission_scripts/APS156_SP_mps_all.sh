#!/bin/bash

#paths
DATE=1201
MSA1=/scratch/aps376/recombo/APS156_SP_Archive/SP_MASTER_OUT/MSA_SP_MASTER
MSA2DIR=/scratch/aps376/recombo/APS156_SP_Archive/500chunks
OUTDIR=/scratch/aps376/recombo/APS156_SP_Archive/${DATE}_SP_mps_dists
SUBMITDIR=/scratch/aps376/recombo/APS156_sp_mps/${DATE}_submissions
SLURMDIR=/scratch/aps376/recombo/APS156_sp_mps/${DATE}_slurm

mkdir -p $OUTDIR
mkdir -p $SUBMITDIR
mkdir -p $SLURMDIR

for i in 0
#for i in {0..23}
do
  echo "submitting part $i"
  jobfile=${SUBMITDIR}/mps_${i}.sh
  MSA2=${MSA2DIR}/MSA_chunk${i}
  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=SP_MPS_$i
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_${DATE}_SP_MPS.out

module load git/gnu/2.16.2
module load go/1.13.6 #try go/1.13.6
module load python3/intel/3.6.3 ##do 3.7.3!
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
#module load sra-tools/intel/2.9.6 #try 2.9.6
module load sra-tools/2.10.5
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

#mcorr
export PATH=\$PATH:\$HOME/go/bin:\$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=\$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=\$PATH:~/opt/ReferenceAlignmentGenerator

cd $OUTDIR
echo \"let's rock\"
mcorr-pair-sync $MSA1 $MSA2 ${OUTDIR}/SP_MPS_${i}_XMFA_OUT.csv  --max-corr-length=3" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done
