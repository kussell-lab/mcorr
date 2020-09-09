#!/bin/bash
##first write for one replicate
##then write for multiple

##make job directory
JOBDIR=$SCRATCH/APS137_ngs_sra

mkdir -p $JOBDIR

for line in 'ngonorrhoeae'
do
  jobfile=$JOBDIR/${line}_sra.sh
  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=${line}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=$JOBDIR/${line}_slurm_%j.out

module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.7.3 ##do 3.7.3!
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
module load sra-tools/2.10.5
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

line=ngonorrhoeae

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path

#mcorr
export PATH=\$PATH:\$HOME/go/bin:\$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=\$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=\$PATH:~/opt/ReferenceAlignmentGenerator

python3 $HOME/ReferenceAlignmentGenerator/FetchSRA_APS2.py $SCRATCH/NGS_SRR_Acc_List.txt $SCRATCH/$line" > $jobfile
  sbatch "$jobfile"
  sleep 1
done
