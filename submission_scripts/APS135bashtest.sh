#!/bin/bash

#SRC is the directory your SRA accession files and genome files are within.
Archive='/Users/asherpreskasteinberg/Desktop/code/recombo/APS135_salmonella'

##MSA stands for multi sequence alignment in the below
SLURM_ARRAY_TASK_ID=1
pipe_dream () {
  #the '$1' command tells it to grab the argument of pipe_dream
	SERO=$1
	#read SERO
	IFS='_'
	read -a strarr <<< "$1"
	echo "main sero: ${strarr[0]}"
#	echo files_$SERO
#	mkdir ${Archive}/${SERO}_OUT
#	cd ${Archive}/${SERO}_OUT

        #ReferenceAlignmentGenerate ${SRC}/SRA_files/sra_accession_$sero ${WRKD}/$sero ${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.fna ${SRC}/Reference/GCF_000006945.2_ASM694v2_genomic.gff ${Archive}/${SERO}_OUT/MSA_$SERO
        #CollectGeneAlignments ${SRC}/SRA_files/sra_accession_$sero
}

echo "SLURM_ARRAY_TASK_ID: " $SLURM_ARRAY_TASK_ID
#this line grabs the line from the list of serotypes in sero_list
#for example when SLURM_ARRAY_TASK_ID=1, you grab the first line
#sero_list=${Archive}/sero_list
sero=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sero_list)
echo "Running Serotype: " $sero
pipe_dream $sero
