#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
DATE=0115
JOBDIR=$SCRATCH/recombo/APS160_sp_allseqs
ARCHIVE=$SCRATCH/recombo/APS160_SP_Archive
#OUTDIR=$SCRATCH/recombo/APS160_SP_Archive/corethreshold95
SUBMITDIR=${JOBDIR}/${DATE}_submissions
SLURMDIR=${JOBDIR}/${DATE}_slurm
mkdir -p ${SUBMITDIR}
mkdir -p ${SLURMDIR}

#mkdir -p $OUTDIR
##list of clusters

##will get to 10 ... in due time
thresholds=(99 95 90)
for line in "${thresholds[@]}"
#for line in 8
do
  for gt in 'CORE' 'FLEX'
  do
    echo "submitting $gt for threshold ${line}"
    jobfile=${SUBMITDIR}/${gt}_threshold_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --mem=30GB
#SBATCH --job-name=${gt}_$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_${gt}_threshold_$line.out

module load go/1.15.2
module load python/intel/3.8.6

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

#mcorr
export PATH=\$PATH:\$HOME/go/bin:\$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=\$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=\$PATH:~/opt/ReferenceAlignmentGenerator

cd $OUTDIR

echo \"let's rock\"
mcorr-xmfa $ARCHIVE/corethreshold${line}/MSA_${gt} $ARCHIVE/corethreshold${line}/corethreshold${line}_${gt}_ALLSEQS_XMFA_OUT &&
mcorr-fit $ARCHIVE/corethreshold${line}/corethreshold${line}_${gt}_ALLSEQS_XMFA_OUT.csv $ARCHIVE/corethreshold${line}/corethreshold${line}_${gt}_ALLSEQS_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
    done
done

