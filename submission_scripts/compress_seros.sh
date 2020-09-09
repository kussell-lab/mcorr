#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS132_compress
SAVEDIR=$ARCHIVE/salmonella/Typhimurium
SCRATCH=$SCRATCH/Typhimurium

mkdir -p $JOBDIR
mkdir -p $SAVEDIR
##make it pause a second in between submissions so we don't break the cluster
for y in 'SRR1957724-SRR1959727' 'SRR1969400-500' 'SRR1969700-SRR1970000' 'SRR1963318-SRR1966998' 'SRR1969500-600' 'SRR1970000-SRR1970327' 'SRR1967022-SRR1969387' 'SRR1969600-700' 'SRR3048518-3049975'
do
  jobfile=$JOBDIR/compress_$y.sh

  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=$y
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/compress_${y}_run-%j.out

module load rclone/1.38

echo \"let's rock\"
tar zcvf $SCRATCH/$y.tar.gz $SCRATCH/$y
rclone copy $y.tar.gz googledrive:hpc/APS132_salmonella_coreflex/Typhimurium
mv $SCRATCH/$y.tar.gz $SAVEDIR/$y.tar.gz" > $jobfile
  sbatch "$jobfile"
  echo "I'm taking a 1 second break"
  sleep 1 #pause the script for a second so we don't break the cluster with our magic
done

