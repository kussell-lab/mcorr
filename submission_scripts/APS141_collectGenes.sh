#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS141_1008_staph_MSA
OUTDIR=$SCRATCH/APS141_Archive
SRC=$HOME/APS136_staph
FASTA=$SRC/Reference/GCF_900004855.1_BB155_genomic.fna
GFF=$SRC/Reference/GCF_900004855.1_BB155_genomic.gff
WRKDIR=$SCRATCH/APS141_saureus
SRA=${WRKDIR}/pileuplist

mkdir -p $JOBDIR
mkdir -p $OUTDIR
##list of clades
for line in "APS141_MSA"
do
  echo "submitting $line"
  jobfile=$JOBDIR/$line.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=4:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_$line.out

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

mkdir $OUTDIR
cd $OUTDIR

echo \"let's rock\"
CollectGeneAlignments ${SRA} $GFF $WRKDIR $OUTDIR/APS141_saureus_MASTER_MSA --appendix \".pileup.fasta\" --progress" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

