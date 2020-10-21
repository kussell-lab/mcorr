
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



##for some reason gene 1826 appears twice, with two different sequences ...
##let's see if any others do

###delete the funny gene
##store line numbers for future deletion
#funnygenes = []

##remove the gene lines

# MSA_file = open(MSA, 'r')
# lines = MSA_file.readlines()
# MSA_file.close()


##re-write the file minus this line
new_file = open(MSA, 'w+')

lastline = lines[0]
for line in lines:
    ##if the line is a gene header, grab the gene name
    ##if it's not the bad gene, write it
    if line.startswith(">"):
        genename = line.split(' ')[0]
        if genename != '>NC_003197.2|cds-YP_009325922.1':
             new_file.write(line)
    ##check if the last line is the bad gene name
    ##if it is, we don't write it
    elif lastline.startswith(">") and genename != '>NC_003197.2|cds-YP_009325922.1':
        new_file.write(line)
    ##check that we don't get two equals in a row
    elif lastline.startswith('=') and line.startswith('=') is False:
        new_file.write(line)
    elif line.startswith('='):
        new_file.write(line)
    lastline = line

##close the file, no moar writing
new_file.close()



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
    #print(genename)
    for ln in master_file:
        ##add a count for the first line
        genelen = genelen + 1
        ##when a new gene begins
        if lastline.startswith(">"):
            thisgene = lastline.split(' ')[0]
            if thisgene != genename:
                print(thisgene)
                genename = thisgene
        if lastline.startswith("="):
            print(lastline)
            print(ln)
        if lastline.startswith("=") and ln.startswith("="):
            print('double equals')
            print(str(i))
            #if the length of the previous gene was
            ##greater than the total number of strains + 1
            ##print the genename
            if genelen > 4710*2+1:
                print('long gene:')
                print(str(i))
                print(genename)
            ##reset the genename for this gene
            genename = ln.split(' ')[0]
            ##reset the gene length count
            genelen = 0
        lastline = ln
        i = i + 1
##check if anything is missing a header
with open(MSA, 'r') as master_file:
    i = 0
    num = 1
    lastline = 'blargh'
    for ln in master_file:
        if lastline.startswith("="):
            num = num + 1
            if not ln.startswith(">"):
                print('missing header: ')
                print(str(i))
                print(ln)
                funnygenes.append(i)

##other crap

# MSA_file = open(MSA, 'r')
# lines = MSA_file.readlines()
# MSA_file.close()
#
#
# ##re-write the file minus this line
# new_file = open(MSA, 'w+')
#
# firstline = '>NC_003197.2|cds-NP_461863.1 3081922-3082482 SRR1968068'
# ##set the name of the current gene
# currentgene = firstline.split(' ')[0]
# for line in lines:
#     ##if the line is a gene header, grab the gene name
#     ##if it's not the bad gene, write it
#     if line.startswith(">"):
#         gene = line.split(' ')[0]
#         ##if the gene header isn't equal to the current one
#         ## put in an equals sign
#         if gene != currentgene:
#             ##put in equals sign
#             new_file.write("=\n")
#             #write the line
#             new_file.write(line)
#             #reset current gene to the new gene
#             currentgene = gene
#
# ##close the file, no moar writing
# new_file.close()