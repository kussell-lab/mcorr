##make a list of pileups to be completed
import os
import numpy as np

dir = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS150_spneumoniae/'
outdir = dir+'tbc/'
complete = set()
piles = open(dir + 'APS150_completepiles', 'r')
for _, pile in enumerate(piles):
    pile = str.rstrip(pile)
    complete.add(pile)

all_fastqs = set()

dwnlds = open(dir + 'APS150_downloaded_sra_list', 'r')
for _, sra in enumerate(dwnlds):
    sra = str.rstrip(sra)
    all_fastqs.add(sra)

tbc = all_fastqs.difference(complete)

tbc = list(tbc)

split_tbc = np.array_split(tbc, 10)

if not os.path.exists(outdir):
    os.makedirs(outdir)
for i in np.arange(0, len(split_tbc)):
    tbc_i = split_tbc[i]
    with open(outdir+'pilesup_TBC_'+str(i), 'w+') as pilesupTBC:
        for sra in tbc_i:
            pilesupTBC.write(sra+'\n')

print(len(tbc))