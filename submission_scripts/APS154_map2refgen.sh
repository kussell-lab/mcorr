#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
DATE=1122
JOBDIR=$BEEGFS/APS154_map2refgen
WRKDIR=$SCRATCH/recombo/saureus
OUTDIR=$BEEGFS/saureus
SRC=$HOME/APS136_staph
FASTA=$SRC/Reference/GCF_900004855.1_BB155_genomic.fna
LISTS=$JOBDIR/piles_tbc3
SUBMITDIR=$JOBDIR/${DATE}_submissions
SLURMDIR=${DATE}_slurm

mkdir -p $SUBMITDIR
mkdir -p $SLURMDIR


##will change to 0 to 9 once confirmed that it werks
for line in {0..24}
do
  echo "submitting list ${line}"
  jobfile=$SUBMITDIR/APS154_map2ref_${line}.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=6:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=map2ref_${line}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_map2ref_${line}.out

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
bash  MapRead2Reference.sh ${LISTS}/pilesup_TBC_${line} $WRKDIR $FASTA" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

