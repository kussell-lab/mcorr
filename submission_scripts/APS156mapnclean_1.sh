#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
DATE=1130_scratch_v1
JOBDIR=$SCRATCH/recombo/APS156_mapnclean
WRKDIR=$SCRATCH/recombo/spneumoniae
OUTDIR=$SCRATCH/recombo/spneumoniae
#WRKDIR=$BEEGFS/spneumoniae
#OUTDIR=$BEEGFS/spneumoniae
SRC=/home/aps376/APS150_spneumoniae
FASTA=$SRC/Reference/GCF_000007045.1_ASM704v1_genomic.fna
#LISTS=$JOBDIR/beegfs_piles1
LISTS=$JOBDIR/scratch_piles_tbc4
SUBMITDIR=${DATE}_submissions
SLURMDIR=${DATE}_slurm

mkdir -p $SUBMITDIR
mkdir -p $SLURMDIR


##will change to 0 to 9 once confirmed that it werks
for line in {2..109}
#for line in 2
do
  echo "submitting list ${line}"
  jobfile=$SUBMITDIR/APS156mapnclean_${line}.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --time=0:45:00
#SBATCH --mem=2GB
#SBATCH --job-name=mapnclean_${line}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_mapnclean_${line}.out

module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
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

##set perl language variable; this will give you fewer warnings
export LC_ALL=C

cd $WRKDIR

echo \"let's rock\"
mapnclean ${LISTS}/piles_TBC_${line} $WRKDIR $FASTA" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 2 second break"
  sleep 2 #pause the script for a second so we don't break the cluster with our magic
done

