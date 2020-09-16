#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS138_cjejuni_MSA
WRKDIR=$SCRATCH/cjejuni
OUTDIR=$SCRATCH/APS138_Archive
SRC=$HOME/APS138_cjejuni
FASTA=$SRC/Reference/GCF_000009085.1_ASM908v1_genomic.fna
GFF=$SRC/Reference/GCF_000009085.1_ASM908v1_genomic.gff
SRA=$HOME/APS138_cjejuni/filereport_read_run_PRJEB31119_tsv.txt

mkdir -p $JOBDIR
mkdir -p $OUTDIR
mkdir -p $WRKDIR
##list of clades
for line in 'APS138_cjejuni'
do
  echo "submitting $line"
  jobfile=$JOBDIR/${line}.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=${line}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_${line}.out

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

mkdir $OUTDIR/${line}_OUT
cd $WRKDIR

echo \"let's rock\"
ReferenceAlignmentGenerateENA ${SRA} $WRKDIR $FASTA $GFF $OUTDIR/${line}_OUT/MSA_$line" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

