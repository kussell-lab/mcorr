import glob
import os

pileupList = glob.glob('*.pileup.fasta')
tarList = glob.glob('*.tar.gz')
file = open('APS156_completepiles_tars', 'w+')
i = 0
for pileup in pileupList:
    size = os.stat(pileup).st_size
    if size != 0:
        i = i + 1
        pilestr = pileup.split('.')
        file.write(pilestr[0]+"\n")
for tar in tarList:
    size = os.stat(tar).st_size
    if size != 0:
        i = i + 1
        tarstr = tar.split('.')
        file.write(tarstr[0]+"\n")
file.close()
print(str(i))

import glob
import os

##also remove ERR907282