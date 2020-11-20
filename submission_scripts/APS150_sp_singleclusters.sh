#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/recombo/APS150_sp_singleclusters
OUTDIR=$SCRATCH/recombo/APS150_SP_Archive

#mkdir -p $OUTDIR
##list of clades

##will get to 10 ... in due time
clusters=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 17 20 21 22 23 24 25 27 28
29 30 31 32 34 35 36 37 38 39 40 41 43 44 45 46 47)
for line in "${clusters[@]}"
do
  for gt in 'CORE' 'FLEX'
  do
    echo "submitting $gt for cluster ${line}"
    jobfile=${gt}_cluster_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=${gt}_cluster$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=slurm%j_${gt}_cluster_$line.out

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

cd $OUTDIR

echo \"let's rock\"
mcorr-xmfa $OUTDIR/cluster${line}/MSA_${gt}_cluster${line} $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT &&
mcorr-fit $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT.csv $OUTDIR/cluster${line}/cluster${line}_${gt}_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
    done
done

