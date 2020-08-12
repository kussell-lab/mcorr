#!/bin/bash
# intro_run.sh

#NOTICE!!!
#You want it to run in this shell session, so this script should be run like this:
#source intro_run.sh

#BACKGROUND
#This script installs all of the tools needed to run the mcorr pipeline using 'ReferenceAlignmentGenerator'
#This will only work if you are on the NYU cluster

#Load all the modules that are needed for ReferenceAlignmentGenerator
echo "Loading modules."
module load git/gnu/2.16.2
module load go/1.10.2 #try go/1.13.6
module load python3/intel/3.6.3 ##do 3.7.3!
module load parallel/20171022
module load prokka/1.12
module load muscle/intel/3.8.31
module load sra-tools/intel/2.8.1-2 #try 2.9.6
module load samtools/intel/1.6
module load smalt/intel/0.7.6
alias roary='singularity exec /beegfs/work/public/singularity/roary-20181203.simg roary'

##This is grabbing dependencies for everything.
echo "Grabbing dependencies."
go get -u github.com/apsteinberg/mcorr/cmd/mcorr-xmfa
go get -u github.com/apsteinberg/mcorr/cmd/mcorr-bam
cd $HOME/go/src/github.com/apsteinberg/mcorr/cmd/mcorr-fit
python3 setup.py install --user
cd ~/

pip install --user tqdm biopython
go get -u github.com/cheggaaa/pb
go get -u github.com/mattn/go-sqlite3
go get -u gopkg.in/alecthomas/kingpin.v2
go get -u github.com/kussell-lab/biogo/seq

#Grabbing ReferenceAlignmentGenerator dependencies.
go get -u github.com/kussell-lab/go-misc/cmd/GenomicConsensus
go get -u github.com/kussell-lab/go-misc/cmd/CollectGeneAlignments

#Grabbing the ReferenceAlignmentGenerator and AssemblyAlignmentGenerator
git clone https://github.com/kussell-lab/AssemblyAlignmentGenerator ~/opt/AssemblyAlignmentGenerator
git clone https://github.com/kussell-lab/ReferenceAlignmentGenerator.git ~/opt/ReferenceAlignmentGenerator

##Making the AssemblyAlignmentGenerator and ReferenceAlignmentGenerator in path
echo "Making everything in path."
#mcorr
export PATH=$PATH:$HOME/go/bin:$HOME/.local/bin

#ReferenceAlignmentGenerator
export PATH=$PATH:~/opt/AssemblyAlignmentGenerator/
export PATH=$PATH:~/opt/ReferenceAlignmentGenerator

echo "All done."
