#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
DATE=0122_single_90
JOBDIR=$SCRATCH/recombo/APS160.5_lmfit
OUTDIR=$SCRATCH/recombo/APS160_SP_Archive/corethreshold90
SUBMITDIR=${JOBDIR}/${DATE}_submissions
SLURMDIR=${JOBDIR}/${DATE}_slurm
mkdir -p ${SUBMITDIR}
mkdir -p ${SLURMDIR}

#mkdir -p $OUTDIR
##list of clusters

##will get to 10 ... in due time
clusters=(8 9 27 75 83 85 89 94 99 106 110 111 112 133 136 140 141 142 149
152 155 158 159 161 162 163 165 169 170 171 173 174 178 180 183 191 196
198 199 213 216 218 220 221 228 229)
for line in "${clusters[@]}"
#for line in 8
do
  for gt in 'CORE' 'FLEX'
  do
    echo "submitting $gt for cluster ${line}"
    jobfile=${SUBMITDIR}/${gt}_cluster_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:20:00
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
fit-stats $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT.csv $OUTDIR/cluster${line}/cluster${line}_${gt}_0122_newparams1" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
    done
done

