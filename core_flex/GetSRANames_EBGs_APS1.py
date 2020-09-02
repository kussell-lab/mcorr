import pandas as pd
import numpy as np
import os
import errno

#original scripts by Spencer, with modifications by Asher

#Background
#This script takes the supplementary table from the Ashton 2016 paper
#and produces a series of SRA accession list files that can be used as ReferenceAlignmentGenerator inputs
#These files are assumed to be generated for the rest of the pipeline, and placed in a folder called 'SRA_files'.

#Specify where the file is that specififies the SRA reference (Ashton 2016 Paper)
supp_path = '~/Desktop/code/peerj_salmonella/SupplementaryTable1.1.xlsx'
#path for outputs
output_path='/Users/asherpreskasteinberg/Desktop/code/recombo/APS135_salmonella/'
#generate sra file folder
sra='SRA_files/'
srapath = os.path.join(output_path, sra)

os.makedirs(srapath, exist_ok=True)


#The prefix that each SRA file with have
out_name_prefix = 'sra_accession'
#The number of serotypes to get
sample_number = 11

#Read the file
Samples_Database = pd.read_excel(supp_path)

#Make column names sane
Samples_Database.columns = [c.replace('-',' ').replace(' ', '_') for c in Samples_Database.columns]

#Pull the Subspecies I rows
Samples_Database_SEI = Samples_Database[Samples_Database.K_mer_identification.str.contains('Salmonella enterica subsp I enterica')]
totcount=len(Samples_Database_SEI)
print(totcount)
#Grab the top represented serotypes
sero_counts = pd.value_counts(Samples_Database_SEI['MLST_derived_serotype'])
##divide the MLST derived serotypes into their eBGs!
#sero_counts = pd.value_counts(Samples_Database_SEI['E_Burst_Group'])
#sero_counts = pd.value_counts(Samples_Database_SEI['Phenotypic_serotype'])
top_sero = sero_counts[0:sample_number].index
EBG_list = []
i = 0
for filter_serotype in top_sero:

    #Grab one serotype
    Samples_Database_SEI_Sub = Samples_Database_SEI[Samples_Database_SEI['MLST_derived_serotype'] == filter_serotype]

    for E_Burst_Group in Samples_Database_SEI_Sub.groupby('E_Burst_Group'):
        #print(E_Burst_Group[0])
        #print(E_Burst_Group[0])
        #Get a row of accession codes
        #out_df = pd.DataFrame([sra for sra in Samples_Database_SEI_Sub.iloc[:,-1] if sra != 'None'])
        out_df = pd.DataFrame([sra for sra in E_Burst_Group[1].iloc[:,-1] if sra != 'None'])
        if len(out_df) < 10:
            continue
        EBGname=filter_serotype.replace('/','-').replace(' ','_')+'-EBG'+str(E_Burst_Group[0])
        out_df.to_csv(output_path + sra + out_name_prefix + '_' + EBGname, sep='\t',index=False,header=False)
        num = len(out_df)
        print('Number of strains in EBG %s of %s is %s' %(E_Burst_Group[0], filter_serotype, num))
        EBG_list.append(EBGname)
        i = i+1
#sra_df = pd.DataFrame([sero.replace('/','-').replace(' ','_') for sero in top_sero])
#sra_df.to_csv(output_path+'sero_list', sep='\t',index=False,header=False)
#EBG_list.to_csv(output_path+'sero_list', sep='\t',index=False,header=False)
with open(output_path+'sero_list','w+') as fh:
    for name in EBG_list:
        fh.write(str(name)+"\n")
print('total count is %s' %i)