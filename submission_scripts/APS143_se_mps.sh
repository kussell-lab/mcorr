#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
SCRATCH=/scratch/aps376
SLURMOUT=${SCRATCH}/recombo/APS143_se_mps/1023_slurm
SUBMISSION=${SCRATCH}/recombo/APS143_se_mps/1023_submissions
OUTDIR=${SCRATCH}/recombo/APS143_1008_senterica_Archive/dists
MSA=${SCRATCH}/recombo/APS143_1008_senterica_Archive/MSA_Master_Sorted
names=${SCRATCH}/recombo/APS143_se_pairs/1023_pairs/1023_se

mkdir -p $SCRATCH
mkdir -p $OUTDIR
mkdir -p $SLURMOUT
mkdir -p $SUBMISSION
##eventually do 0 to 79

for line in {0..1}
do
  echo "submitting job #${line}"
  #add $JOBDIR
  jobfile=${SUBMISSION}/APS143_SE_MPS_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=se_$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMOUT}/slurm%j_se_mps_$line.out

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
mcorr-pair-specific $MSA $OUTDIR/HP_mps_${line}_XMFA_OUT.csv --pair-list=${names}_${line}.csv &&
mcorr-fit $OUTDIR/SE_MPS_${line}_XMFA_OUT.csv $OUTDIR/SE_MPS_${line}_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

