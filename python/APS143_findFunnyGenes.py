
'''
Goal:
see if there is a gene that has not been properly
annotated in the APS143 s enterica  MSA file
'''


from tqdm import tqdm
import pandas as pd

MSA = "/scratch/aps376/recombo/APS143_1008_senterica_Archive/MSA_Master_Sorted"

##read the file
j = 0
while j == 0:
    funnygenes = []
    with open(MSA, 'r') as master_file:
        i = 0
        lastline = 'blargh'
        for ln in master_file:
            if lastline.startswith("="):
                if not ln.startswith(">"):
                    print(str(i))
                    print(ln)
                    funnygenes.append(i)
            lastline = ln
            i = i + 1
    if len(funnygenes) == 0:
        j = 1
        break

    MSA_file = open(MSA, 'r')
    lines = MSA_file.readlines()
    MSA_file.close()

    ##delete the line in question
    for gene in funnygenes:
        del lines[gene]

    ##re-write the file minus this line
    new_file = open(MSA, 'w+')

    for line in lines:
        new_file.write(line)

    new_file.close()