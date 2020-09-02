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
#specify where the file lives which has a list of MLSTs and clonal complexes
cc_path = '~/Desktop/code/recombo/APS136_staph/MLST_profiles.csv'
#specify the path to the SRA numbers so we can go from ERS to SRA
sranumspath = '~/Desktop/code/recombo/APS136_staph/SraRunTable.csv'
#path for outputs
output_path='/Users/asherpreskasteinberg/Desktop/code/recombo/APS136_staph/'
#generate sra file folder
sra='SRA_files/'
srapath = os.path.join(output_path, sra)

os.makedirs(srapath, exist_ok=True)


#The prefix that each SRA file with have
out_name_prefix = 'sra_accession'
#The number of serotypes to get
sample_number = 7

#Read the file
Samples_Database = pd.read_csv(supp_path)

cc_file = pd.read_csv(cc_path, sep='\t')

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

this = cc_file[cc_file['ST'] == Samples_Database['MLST'][11]]
cc = this.iloc[0]['clonal_complex']
ers = Samples_Database['ERS_Number'][11]
row = sraruntable[sraruntable['SRA_accession'].str.match(ers)]
row = row[row.Bases == row.Bases.max()]
err = row.iloc[0]['Run']
data = [{'clonal_complex': cc, 'MLST': Samples_Database['MLST'][11], 'ERS_Number':ers, 'ERR_Number':err}]
staph = pd.DataFrame(data)

for i in np.arange(12, len(Samples_Database)):
    ers = Samples_Database['ERS_Number'][i]
    if np.isnan(Samples_Database['MLST'][i]) == True:
        continue
    if ers == '#N/A':
        continue
    if ers not in sraruntable.values:
        continue
    this = cc_file[cc_file['ST'] == Samples_Database['MLST'][i]]
    cc = this.iloc[0]['clonal_complex']
    row = sraruntable[sraruntable['SRA_accession'].str.match(ers)]
    row = row[row.Bases == row.Bases.max()]
    err = row.iloc[0]['Run']
    data = [{'clonal_complex': cc, 'MLST': Samples_Database['MLST'][i], 'ERS_Number':ers, 'ERR_Number':err}]
    staph = staph.append(data)
staph = staph.dropna()
staph = staph.reset_index(drop=True)
staph.to_csv(output_path + 'staph_df.csv')



cc_counts = pd.value_counts(staph['clonal_complex'])


#Pull the Subspecies I rows
#Samples_Database_SEI = Samples_Database[Samples_Database.K_mer_identification.str.contains('Salmonella enterica subsp I enterica')]

#Grab the top represented clonal complexes
#top_cc = cc_counts[0:sample_number].index
cc_names = []

for clonal_complex in staph.groupby('clonal_complex'):
    if len(clonal_complex[1]) < 10:
        continue
    cc = clonal_complex[1]['ERR_Number']
    cc.to_csv(output_path + sra + out_name_prefix + '_' + clonal_complex[0], sep='\t',index=False,header=False)
    cc_names.append(clonal_complex[0])

with open(output_path+'cc_list', 'w') as f:
    i = 1
    for item in cc_names:
            f.write(item + '\n')


# for filter_serotype in top_sero:
#
#     #Grab one serotype
#     Samples_Database_SEI_Sub = Samples_Database_SEI[Samples_Database_SEI['MLST_derived_serotype'] == filter_serotype]
#
#     #Get a row of accession codes
#     out_df = pd.DataFrame([sra for sra in Samples_Database_SEI_Sub.iloc[:,-1] if sra != 'None'])
#     out_df.to_csv(output_path + sra + out_name_prefix + '_' + filter_serotype.replace('/','-').replace(' ','_'), sep='\t',index=False,header=False)
#
# sra_df = pd.DataFrame([sero.replace('/','-').replace(' ','_') for sero in top_sero])
# sra_df.to_csv(output_path+'sero_list', sep='\t',index=False,header=False)
