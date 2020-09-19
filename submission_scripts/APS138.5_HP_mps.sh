#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
SCRATCH=/scratch/aps376
JOBDIR=${SCRATCH}/APS138.5_HP_mps
OUTDIR=${SCRATCH}/APS138.5_Archive
SRC=$HOME/APS138.5
MSA=${SCRATCH}/Helicobacter_pylori_global_population/BIGSdb_gene-by-gene_alignment.xmfa
names=${SCRATCH}/Helicobacter_pylori_global_population/strain_lists/0918_hpylori

mkdir -p $SCRATCH
mkdir -p $JOBDIR
mkdir -p $OUTDIR
##eventually do 0 to 79

for line in {0..2}
do
  echo "submitting job #${line}"
  #add $JOBDIR
  jobfile=$JOBDIR/APS138.5_HP_mps_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=HP_$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_HP_mps_$line.out

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
mcorr-fit $OUTDIR/HP_mps_${line}_XMFA_OUT.csv $OUTDIR/HP_mps_${line}_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

