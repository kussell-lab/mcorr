#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
DATE=0113
JOBDIR=$SCRATCH/recombo/APS159_ng_singleclusters
OUTDIR=$SCRATCH/recombo/APS159_NG_Archive
SUBMITDIR=${JOBDIR}/${DATE}_submissions
SLURMDIR=${JOBDIR}/${DATE}_slurm
mkdir -p ${SUBMITDIR}
mkdir -p ${SLURMDIR}

#mkdir -p $OUTDIR
##list of clusters

##will get to 10 ... in due time
clusters=(12 36 40 43 57 60)
for line in "${clusters[@]}"
#for line in 3
do
  for gt in 'CORE'
  do
    echo "submitting $gt for cluster ${line}"
    jobfile=${SUBMITDIR}/${gt}_cluster_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=${gt}_cluster$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_${gt}_cluster_$line.out

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
mcorr-xmfa $OUTDIR/cluster${line}/MSA_${gt}_cluster${line} $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT &&
mcorr-fit $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT.csv $OUTDIR/cluster${line}/cluster${line}_${gt}_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
    done
done

