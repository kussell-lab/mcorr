import glob
import os

tarList = glob.glob('*_1.fastq.gz')
file = open('APS156_scratch_zips', 'w+')
i = 0
for tar in tarList:
    # size = os.stat(tar).st_size
        tarstr = tar.split('_1.')
        file.write(tarstr[0]+"\n")
        i = i+1
file.close()
print(str(i))