#!/bin/bash

#paths
DATE=0309
ARCHIVE=/scratch/aps376/recombo/APS170_SA_Archive
MSA1=${ARCHIVE}/MSA_SA_MASTER_GAPFILTERED
MSA2DIR=${ARCHIVE}/1000chunks
OUTDIR=${ARCHIVE}/${DATE}_SA_mps_dists
SUBMITDIR=/scratch/aps376/recombo/APS170clusterseqs/${DATE}_submissions
SLURMDIR=/scratch/aps376/recombo/APS170clusterseqs/${DATE}_slurm

mkdir -p $OUTDIR
mkdir -p $SUBMITDIR
mkdir -p $SLURMDIR

#for i in 0
for i in {0..1004}
do
  echo "submitting part $i"
  jobfile=${SUBMITDIR}/mps_${i}.sh
  MSA2=${MSA2DIR}/MSA_chunk${i}
  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=12:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=SA_MPS_$i
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_${DATE}_SA_MPS.out

module load go/1.15.7

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

#mcorr
export PATH=\$PATH:\$HOME/go/bin:\$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=\$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=\$PATH:~/opt/ReferenceAlignmentGenerator

cd $OUTDIR
echo \"let's rock\"
mcorr-pair-sync $MSA1 $MSA2 ${OUTDIR}/SA_MPS_${i}_XMFA_OUT.csv  --max-corr-length=3" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done
