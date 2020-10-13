import math

import pandas as pd
import numpy as np
import os
import errno

#original scripts by Spencer, with modifications by Asher

#Background
#This script takes the supplementary table from the Ashton 2016 paper
#and produces a series of SRA accession list files that can be used as ReferenceAlignmentGenerator inputs
#These files are assumed to be generated for the rest of the pipeline, and placed in a folder called 'SRA_files'.

#Specify where the file is that specififies the ERS numbers
supp_path = '~/Desktop/code/recombo/APS136_staph/microreact-project-EkUvg9uY-data.csv'

#specify the path to the SRA numbers so we can go from ERS to SRA
sranumspath = '~/Desktop/code/recombo/APS136_staph/SraRunTable.csv'
#path for outputs
output_path='/Users/asherpreskasteinberg/Desktop/code/recombo/APS141/'


#The prefix that each SRA file with have
out_name_prefix = 'sra_accession'

#Read the file
Samples_Database = pd.read_csv(supp_path)

##read the SRA accession list
sraruntable = pd.read_csv(sranumspath)
sraruntable = sraruntable[sraruntable['Organism'] == 'Staphylococcus aureus']
# ers = staph['ERS_Number']
# run_nums = []
# for i in np.arange(0, len(ers)+1):
#     print(i)
#     row = sraruntable[sraruntable['SRA_accession'].str.match(ers[i])]
#     row = row[row.Bases == row.Bases.max()]
#     err = row['Run']
#     run_nums.append(row.iloc[0]['Run'])

#Make column names sane
Samples_Database.columns = [c.replace('-',' ').replace(' ', '_') for c in Samples_Database.columns]

print(len(Samples_Database))

##match the ers number to the sra accession number
ers = Samples_Database['ERS_Number'][5]
row = sraruntable[sraruntable['SRA_accession'].str.match(ers)]
row = row[row.Bases == row.Bases.max()]
err = row.iloc[0]['Run']
data = [{'ERS_Number':ers, 'ERR_Number':err}]
staph = pd.DataFrame(data)

for i in np.arange(6, len(Samples_Database)):
    ers = Samples_Database['ERS_Number'][i]
    if ers == '#N/A':
        continue
    if ers not in sraruntable.values:
        continue
    row = sraruntable[sraruntable['SRA_accession'].str.match(ers)]
    row = row[row.Bases == row.Bases.max()]
    err = row.iloc[0]['Run']
    data = [{'ERS_Number':ers, 'ERR_Number':err}]
    staph = staph.append(data)
staph = staph.dropna()
staph = staph.reset_index(drop=True)
staph.to_csv(output_path + 'staph_df.csv')
staph['ERR_Number'].to_csv(output_path+'sra_accession_list',
                           index = False, header = False)
print(len(staph))