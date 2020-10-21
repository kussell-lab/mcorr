
'''
Goal:
see if there is a gene that has not been properly
annotated in the APS143 s enterica  MSA file
this time looking at where mcorr-pair broke from the slurm output file
'''


from tqdm import tqdm
import pandas as pd

MSA = "/scratch/aps376/recombo/APS143_1008_senterica_Archive/MSA_Master_Sorted"

##read the file

funnygenes = []

with open(MSA, 'r') as master_file:
    i = 0
    num = 1
    lastline = 'blargh'
    for ln in master_file:
        if lastline.startswith("="):
            num = num + 1
            if not ln.startswith(">"):
                print(str(i))
                print(ln)
                funnygenes.append(i)
        ##this is the line that mcorr-pair stopped at
#        if num == 1826:
#            print(ln)
        ##this is the line after
#        if num == 1827:
#            print(ln)
#        lastline = ln
 #       i = i + 1

##for some reason gene 1826 appears twice, with two different sequences ...
##let's see if any others do

###delete the funny gene
##store line numbers for future deletion
#funnygenes = []

# with open(MSA, 'r') as master_file:
#     ##line number count
#     i = 0
#     for ln in master_file:
#         if ln.startswith(">"):
#             genename = ln.split(' ')[0]
#             if genename == '>NC_003197.2|cds-YP_009325922.1':
#                 ##get that extra equals sign, which may screw things up
#                 if lastline.startswith("="):
#                     funnygenes.append(i-1)
#                 #grab the gene name line, and the next line which is the gene
#                 funnygenes.append(i)
#                 funnygenes.append(i+1)
#             ##get our other culprit
#             elif genename == '>NC_003197.2|cds-NP_459150.1':
#                 if lastline.startswith("="):
#                     funnygenes.append(i-1)
#                 funnygenes.append(i)
#                 funnygenes.append(i+1)
#         lastline = ln
#         i = i + 1

##remove the gene lines

# MSA_file = open(MSA, 'r')
# lines = MSA_file.readlines()
# MSA_file.close()
#
# for gene in funnygenes:
#     print(str(gene))
#     del lines[gene]
#
# ##re-write the file minus this line
# new_file = open(MSA, 'w+')
#
# for line in lines:
#     new_file.write(line)
#
# new_file.close()

## double check we gucci now

with open(MSA, 'r') as master_file:
    #line number count
    i = 0
    lastline = 'blargh'
    ##count the length of the genes
    genelen = 0
    ##scan line by line through the file
    ## set the first gene name
    firstline = master_file.readline()
    genename = firstline.split(' ')[0]
    for ln in master_file:
        ##add a count for the first line
        genelen = genelen + 1
        ##when a new gene begins
        if lastline.startswith("="):
            #if the length of the previous gene was
            ##greater than the total number of strains + 1
            ##print the genename
            if genelen > 4710*2+1:
                print('long gene')
                print(str(i))
                print(genename)
            ##reset the genename for this gene
            genename = ln.split(' ')[0]
            ##reset the gene length count
            genelen = 0
        lastline = ln
        i = i + 1

#print(funnygenes)


