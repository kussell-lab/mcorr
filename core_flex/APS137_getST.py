#!/usr/bin/env python
import os
import pandas as pd
import glob
import numpy as np

#define directories
#directory with the MLST csv files
outdir = "/Users/asherpreskasteinberg/Desktop/code/recombo/APS137_ngs/"
ngs = '/Volumes/GoogleDrive/My Drive/hpc/APS137_MLST_output/'
pubmlstpath = "/Users/asherpreskasteinberg/Desktop/code/recombo/APS137_ngs/pubmlst_neisseira.csv"

#compile a dataframe with all the ST profile outputs from stringMLST
os.chdir(ngs)
dat = pd.read_csv('SRR3349835_MLST', sep='\t')

for file in glob.glob("*MLST"):
    #print(file)
    if file == 'SRR5007203_MLST':
        continue
    dat1 = pd.read_csv(file, sep='\t')
    if len(dat1) > 1:
        continue
    dat = dat.append(dat1)
print(len(dat))
dat = dat.reset_index(drop=True)

##put the pubmlst clonal complex list into a dataframe
pubmlst = pd.read_csv(pubmlstpath, sep='\t')

##make a np array for storing the clonal complexes for each dataframe entry
cc = []
##match MLST profile to a clonal complex
for i in np.arange(0, len(dat)):
    #print(i)
    if dat['ST'][i] == 0:
        cc.append('nan')
    else:
        cc_row = pubmlst[pubmlst['ST'] == dat['ST'][i]]
        one_cc = cc_row.iloc[0]['clonal_complex']
        cc.append(one_cc)

# Convert the dataframe to an XlsxWriter Excel object.
print(len(cc))

##save the MLST profiles and clonal complexes
cc_dat = pd.DataFrame(cc, columns=['clonal_complex'])
cc_dat.to_csv(outdir+'APS137_clonal_complex.csv')
dat['clonal_complex'] = cc_dat['clonal_complex']
dat.to_excel(outdir+"APS137_MLST_profiles.xlsx")
dat.to_csv(outdir+'APS137_MLST_profiles.csv')

##count how many strains have the same MLST
MLST_counts = pd.value_counts(dat['ST'])
##first column is the ST, second is counts
MLST_counts.to_csv(outdir+'ST_counts.csv', header = ['count'])

##how many STs include greater than 10 strains?
big_counts = MLST_counts[MLST_counts >= 10]
##put the strains we're including into a .csv file
big_counts.to_csv(outdir+'large_ST_counts.csv', header = ['count'])

#generate sra file folder
sra='SRA_files/'
srapath = os.path.join(outdir, sra)

os.makedirs(srapath, exist_ok=True)


#The prefix that each SRA file with have
out_name_prefix = 'sra_accession'

##group by ST and put into sra accession lists so we can make MSA files
stnames = []
for ST in dat.groupby('ST'):
    if len(ST[1]) < 10:
        continue
    stgroup = ST[1]['Sample']
    stgroup.to_csv(outdir+sra+out_name_prefix+'_ST-'+str(ST[0]), sep='\t',index=False,header=False)
    stnames.append(str(ST[0]))

with open(outdir+'ST_list','w') as f:
    for item in stnames:
        f.write(item+ '\n')