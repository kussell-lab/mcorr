#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS139_ngs_map2ref
WRKDIR=$SCRATCH/ngonorrhoeae
OUTDIR=$SCRATCH/APS139_Archive
SRC=$HOME/APS137_ngs
FASTA=$SRC/Reference/GCF_000020105.1_ASM2010v1_genomic.fna
GFF=$SRC/Reference/GCF_000020105.1_ASM2010v1_genomic.gff
LISTS=$SCRATCH/APS139_ena_lists

mkdir -p $JOBDIR
mkdir -p $OUTDIR

##will change to 0 to 9 once confirmed that it werks
for line in {0..1}
do
  echo "submitting list ${line}"
  jobfile=$JOBDIR/APS139_map2ref_${line}.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=map2ref_${line}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_map2ref_${line}.out

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

cd $WRKDIR

echo \"let's rock\"
bash  MapRead2Reference.sh ${LISTS}/ngs_list_${line} $WRKDIR $FASTA" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

