import glob

pileupList = glob.glob('*.pileup.fasta')
file = open('APS150_completepiles', 'w+')
for pileup in pileupList:
    pilestr = pileup.split('.')
    file.write(pilestr[0]+"\n")
file.close()
print(len(pileupList))

##also remove ERR907282