#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
SCRATCH=/scratch/aps376
JOBDIR=${SCRATCH}/APS139_ngs_mps
OUTDIR=${SCRATCH}/APS139_mps_Archive
MSA=${SCRATCH}/MSA_CORE_Master_Sorted
names=${SCRATCH}/APS139_ngs_strain_lists/0924_ngs

mkdir -p $SCRATCH
mkdir -p $JOBDIR
mkdir -p $OUTDIR
##eventually do 0 to 99

for line in {0..3}
do
  echo "submitting job #${line}"
  #add $JOBDIR
  jobfile=$JOBDIR/APS139_ngs_mps_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=NGS_$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_ngs_mps_$line.out

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
mcorr-pair-specific $MSA $OUTDIR/HP_mps_CORE_${line}_XMFA_OUT.csv --pair-list=${names}_${line}.csv" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

