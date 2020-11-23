import glob
import os

pileupList = glob.glob('*.pileup.fasta')
file = open('APS154_blankpiles', 'w+')
for pileup in pileupList:
    size = os.stat(pileup).st_size
    if size == 0:
        pilestr = pileup.split('.')
        file.write(pilestr[0]+"\n")
file.close()
print(len(pileupList))

##also remove ERR907282