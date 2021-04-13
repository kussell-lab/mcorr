#!/bin/bash
DATE=210331_SA_mcorr
ARCHIVE=/scratch/aps376/recombo/APS170_SA_Archive
OUTDIR=/scratch/aps376/recombo/APS180_SA_Archive
MSA=${ARCHIVE}/MSA_SA_MASTER_GAPFILTERED
list=${ARCHIVE}/strain_list
SUBMITDIR=/scratch/aps376/recombo/APS180genebins/${DATE}_submissions
SLURMDIR=/scratch/aps376/recombo/APS180genebins/${DATE}_slurm

mkdir -p ${SUBMITDIR}
mkdir -p ${SLURMDIR}

bins=(0 20 40 60 80)
for i in "${bins[@]}"
do
  min=$i
  max=$(expr $i + 20)
  echo "submitting ${min}-${max}"
  jobfile=${SUBMITDIR}/bin_${min}-${max}.sh
  echo "#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --mem=4GB
#SBATCH --job-name=bin${min}-${max}
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aps376@nyu.edu
#SBATCH --output=${SLURMDIR}/slurm%j_bin_${min}-${max}.out


##INPUTS

echo \"Loading modules.\"
module load go/1.15.7
module load singularity/3.6.4

##activate venv for mcorr-fit
cd /scratch/aps376/recombo
source venv/bin/activate
export OMP_NUM_THREADS=\$SLURM_CPUS_PER_TASK;

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path
echo \"Making everything in path.\"
#mcorr
export PATH=\$PATH:\$HOME/go/bin:\$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=\$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=\$PATH:~/opt/ReferenceAlignmentGenerator

##set perl language variable; this will give you fewer warnings
export LC_ALL=C

echo \"let's rock\"
cd ${OUTDIR}/bin${min}-${max}
mcorr-xmfa MSA_${min}-${max} MSA_${min}-${max}_xmfa_out &&
mcorr-fit MSA_${min}-${max}_xmfa_out.csv MSA_${min}-${max}_fit_out || true" > $jobfile
    sbatch "$jobfile"
    echo "I'm taking a 1 second break"
    sleep 1 #pause the script for a second so we don't break the cluster with our magic
done


