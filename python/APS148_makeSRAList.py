import glob

pileupList = glob.glob('*_1.fastq')
file = open('sra_list', 'w+')
for pileup in pileupList:
    pilestr = pileup.split('_')
    file.write(pilestr[0]+"\n")
file.close()
print(len(pileupList))