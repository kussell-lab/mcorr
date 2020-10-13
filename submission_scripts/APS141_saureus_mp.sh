#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
SCRATCH=/scratch/aps376
JOBDIR=${SCRATCH}/APS141_saureus_mp
OUTDIR=${SCRATCH}/APS141_1008_Archive
MSA=${SCRATCH}/APS141_1008_Archive/APS141_saureus_MASTER_MSA

mkdir -p $JOBDIR
mkdir -p $OUTDIR
##eventually do 0 to 79

for line in 'saureus_mp_1009' #{0..79}
do
  echo "submitting job #${line}"
  #add $JOBDIR
  jobfile=$JOBDIR/APS139_ngs_mp_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_$line.out

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
mcorr-pair $MSA $OUTDIR/saureus_mp_MASTER_XMFA_OUT.csv  --max-corr-length=3" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

