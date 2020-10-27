#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/recombo/APS148_pa_incluster
OUTDIR=$SCRATCH/recombo/APS148_PA_Archive
MSADIR=${OUTDIR}/cutoff_10pt

mkdir -p $JOBDIR
mkdir -p $OUTDIR
##list of clades

##will get to 10 ... in due time
clusters=(1 2 3 4 5 7 9 11 12 13 14 15 16 19 20 21
22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38
39 43 44 45 47 54 55 56 58 61 62 63 64 65 67 68 69
70 73 74 75 76 77 78 79 81 82 84)
for line in "${clusters[@]}"
do
  for gt in 'CORE' 'FLEX'
  do
    echo "submitting $gt for cluster ${line}"
    jobfile=$JOBDIR/${gt}_cluster_${line}.sh

    echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=${gt}_cluster$line
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/slurm%j_${gt}_cluster_$line.out

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

mkdir -p cluster${line}

echo \"let's rock\"
mcorr-xmfa $MSADIR/MSA_${gt}_cluster${line} $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT &&
mcorr-fit $OUTDIR/cluster${line}/cluster${line}_${gt}_XMFA_OUT.csv $OUTDIR/cluster${line}/cluster${line}_${gt}_FIT_OUT || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
    done
done

