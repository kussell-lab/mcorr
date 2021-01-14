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
clusters=(4 5 6 7 8 9 12 13 14 15 16 17 19 20 21 22 23 26 28 29 31 32 33
34 36 37 39 40 41 43 46 48 49 50 51 52 53 54 57 58 59 60 61 62)
for line in "${clusters[@]}"
#for line in 3
do
  for gt in 'CORE' 'FLEX'
  do
    echo "submitting $gt for cluster ${line}"
    jobfile=${SUBMITDIR}/${gt}_cluster_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
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

