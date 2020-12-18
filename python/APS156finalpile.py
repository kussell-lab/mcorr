import glob
import os

pileupList = glob.glob('*.pileup.fasta')
file = open('APS156_finalpile', 'w+')
i = 0
for pileup in pileupList:
    size = os.stat(pileup).st_size
    if size != 0:
        pilestr = pileup.split('.')
        file.write(pilestr[0]+"\n")
        i = i + 1
file.close()
print(str(i))

##also remove ERR907282