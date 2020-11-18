### makeClusters

Package which takes a multi-sequence alignment (MSA) file generated from the ReferenceAlginmentGenerator package, and
pairwise distances and does the following:\
(1) Puts the pairwise distances into a distance matrix called "distancematrix.npy" and the names of each row in
the order corresponding to the distance matrix in a separate text file called "strains"\
(2) Takes the distance matrix and strain list, and breaks the strains into clusters by cutting the sequence tree\
at a percentile of the pairwise distances of the user's choosing (default is 10). This\
outputs a text file called "stats" with information on the clustering and another text file called\
"cluster_list" which relates strains to cluster number.\
(3) Takes the MSA file of all strains and the "cluster_list" text file and breaks\
the MSA file into MSA files for each cluster and MSA files for both the core and flexible genes\
for each cluster, based on a threshold percentage over which a gene is considered a core gene (default is 90).\
The gene also must be present in all clusters for it to be considered core. In addition, it outputs\
a csv file with the percentage of strains which have the gene in each cluster.\
\
To install the package, take the following steps:\
(1) Pull the makeClusters directory from github or copy it to your computer.\
(2) Use "pip install /path/makeClusters" to install the package (where "path" is the path to the makeClusters folder).\
Alternatively use "pip -e install /path/makeClusters" if you want to edit the package.\
(3) Install the dependencies "mcorr-dm" and "write-cluster-msa" by using the following commands:\
    (a) go get -u github.com/apsteinberg/mcorr/cmd/development/mcorr-dm\
    (b) go get -u github.com/apsteinberg/mcorr/cmd/development/write-cluster-msa

For details on basic usage (e.g., input order) after installation, use the command makeClusters --help