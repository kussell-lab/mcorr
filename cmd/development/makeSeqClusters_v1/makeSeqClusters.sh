#!/bin/bash
# This script generates sequence clusters and divides them into core and flexible genomes
#using the output of mcorr-pair or mcorr-pair-sync
# Inputs:
#   (1) mcorr-pair_csv: the output of pairwise distances from mcorr-pair;
#   (2) working_dir: the working space and output directory
#   (3) cutoff: cutoff percentile of pairwise distances to use for making flat clusters (%)
#   (4) master_msa: master msa file for all strain sequences
#   (5) threshold: threshold percentage above which you're considered a core gene (%)

mcorr-pair_csv=$1
working_dir=$2
cutoff=$3
master_msa=$4
threshold=$5

#mcorr_pair_csv=saureus_mp_MASTER_XMFA_OUT.csv
#working_dir=/Volumes/aps_timemachine/makeSeqClusters_test
#cutoff=10
#master_msa=APS141_saureus_MASTER_MSA
#threshold=90

# Obtain the source file directory.

##start timer
SECONDS=0

cd ${working_dir}

#Step 1: take mcorr-pair output and translate to distance matrix and a list of strains
echo "converting mcorr-pair to matrix ..."
mcorr-dm $mcorr_pair_csv distancematrix

#Step 2: take the distance matrix and list of names
echo "clustering sequences ..."
makeClusters distancematrix.npy strains $cutoff

#Step 3: make sequence cluster MSA files and core and flexible MSA files
echo "writing cluster MSA files ..."
write-cluster-msa $master_msa cluster_list --core-cutoff=$threshold

echo "Done with making clusters, time to boogie"

##timer
if (($SECONDS > 3600)); then
  let "hours=SECONDS/3600"
  let "minutes=(SECONDS%3600)/60"
  let "seconds=(SECONDS%3600)%60"
  echo "Completed in $hours hour(s), $minutes minute(s) and $seconds second(s)"
elif (($SECONDS > 60)); then
  let "minutes=(SECONDS%3600)/60"
  let "seconds=(SECONDS%3600)%60"
  echo "Completed in $minutes minute(s) and $seconds second(s)"
else
  echo "Completed in $SECONDS seconds"
fi
