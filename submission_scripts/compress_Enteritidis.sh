#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS132_compress
SAVEDIR=$ARCHIVE/salmonella/Enteritidis
SCRATCH=$SCRATCH/Enteritidis

mkdir -p $JOBDIR
mkdir -p $SAVEDIR
##make it pause a second in between submissions so we don't break the cluster
for y in 'SRR196_' 'SRR1960082-SRR1965290' 'SRR1966_' 'SRR1967_' 'SRR1968_' 'SRR1969_'
do
  jobfile=$JOBDIR/compress_$y.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=72:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=$y
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/compress_${y}_run-%j.out

module load rclone/1.38

echo \"let's rock\"
tar zcvf $SCRATCH/$y.tar.gz $SCRATCH/$y
rclone copy $y.tar.gz googledrive:hpc/APS132_salmonella_coreflex/Enteritidis" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

