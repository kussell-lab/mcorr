import glob
import os

tarList = glob.glob('*.tar.gz')
file = open('APS156_scratch_complete_tars', 'w+')
i = 0
for tar in tarList:
    # size = os.stat(tar).st_size
    # if size != 0:
        tarstr = tar.split('.')
        file.write(tarstr[0]+"\n")
        i = i+1
file.close()
print(str(i))