import pandas as pd
from itertools import combinations
import numpy as np

#some functions for parsing xmfa files may not be relevant
#depending how you feed in strain names

def parse_description_line(line):
    """
    Parse description line of a FASTA record.
    And return the strain and alignment names.
    """
    line = line.rstrip()[1:] # ignore '>'
    words = line.split(' ')
    #alignment_name = words[len(words)-1]
    alignment_name = words[0]
    strain_name = words[1]
    #strain_name = words[0].split(':')[0]
    return (strain_name, alignment_name)

def read_fasta(filename):
    """
    Read a FASTA or XMFA file.
    And return a list of (strain, alignment, sequence).
    """
    sequence_list = []
    strain_name = None
    alignment_name = None
    sequence = ""
    with open(filename, 'rU') as reader:
        for line in reader:
            line = line.strip()
            if line.startswith('>'): # description line
                if strain_name is not None:
                    sequence_list.append((strain_name, alignment_name, sequence))
                strain_name, alignment_name = parse_description_line(line)
                sequence = ""
            elif line.startswith('=') or len(line) == 0:
                continue
            else:
                sequence = sequence + line
    if strain_name is not None and sequence != "":
        sequence_list.append((strain_name, alignment_name, sequence))
    return sequence_list

"""
begin making pairs
"""
##number of splits you want to make
numsplits = 10
##what you want to name the output csv files
outdir = '/Users/asherpreskasteinberg/Desktop/code/recombo/APS138_cjejuni/'
csvname = 'APS138_enalists/ena_list'

#the original list
oglistpath = "/Users/asherpreskasteinberg/Desktop/code/recombo/APS138_cjejuni/138_cj_ENAlists.txt"
ENA = open(oglistpath, 'r')
content = ENA.read()
strains = content.split("\n")
print(len(strains))
split_list = np.array_split(strains, numsplits)
tot = 0

full_list = pd.DataFrame(data = strains)
full_list.to_csv(outdir+csvname+'_full', index = False, header = False)
for i in np.arange(0, len(split_list)):
    tot = tot + (len(split_list[i]))
    single_split = pd.DataFrame(data = split_list[i])
    single_split.to_csv(outdir+csvname+'_%s' %i, index = False, header = False)

print(tot)

#combolist.to_csv('HP50_test.csv')