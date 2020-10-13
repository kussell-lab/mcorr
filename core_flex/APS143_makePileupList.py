import glob
import numpy as np

pileupList = glob.glob('*.pileup.fasta')

pileuparray = np.array(pileupList)

file = open('pileuplist', 'w+')
for pileup in pileupList:
    pilestr = pileup.split('.')
    file.write(pilestr[0]+"\n")
file.close()
print(len(pileupList))