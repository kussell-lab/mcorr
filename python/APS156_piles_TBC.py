##make a list of pileups to be completed
import os
import numpy as np

dir = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS156_SP_analysis/'
outdir = dir+'scratch_piles_tbc4/'
complete = set()
# piles = open(dir + 'APS156_completepiles_tars', 'r')
# for _, pile in enumerate(piles):
#     pile = str.rstrip(pile)
#     complete.add(pile)

all_fastqs = set()

dwnlds = open(dir + 'APS156_scratch_zips', 'r')
#dwnlds = open(dir + 'APS154_blankpiles', 'r')
for _, sra in enumerate(dwnlds):
    sra = str.rstrip(sra)
    all_fastqs.add(sra)

print("all SRA: %d" % len(all_fastqs))
tbc = all_fastqs.difference(complete)

tbc = list(tbc)

split_tbc = np.array_split(tbc, 110)

if not os.path.exists(outdir):
    os.makedirs(outdir)
for i in np.arange(0, len(split_tbc)):
    tbc_i = split_tbc[i]
    with open(outdir+'piles_TBC_'+str(i), 'w+') as pilesupTBC:
        for sra in tbc_i:
            pilesupTBC.write(sra+'\n')

print("to be completed: %d" % len(tbc))