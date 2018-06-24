#Python Script to compare co-localization across trans-contacts for a gene subset
#Using 10000 permutations to create a background model based on pairwise comparisons
#Set at 500KB resolution
#Takes in background list as an optional parameter of which it selects the genes to randomly permute

import os
import sys
import operator
import aopen
import itertools
import ConfigParser
import random

#New addition - October 16, 2017 - include a background list to select from as an optional parameter

#Method to get the experimental and random permutations for a set of genes, given the gene list and the chromosome/bin information
def getCoLocalizationData(genes,bkg_Dict,gene_Dict,bin_chrs,perms,verbose):
    pairwise = list(itertools.combinations(genes, 2))
    if verbose:
        print pairwise
    colocalizationList = list()
    permutation_Dict = {} # dictionary to store lists of permutation scores
    scoreList = list()
    for i in range(perms): #initialize the dictionary with 1000 lists
        permutation_Dict[i]=list()
    for gene_pair in pairwise:
        gene1 = gene_pair[0]
        gene2 = gene_pair[1]
        bin1 = int(gene_Dict[gene1])
        bin2 = int(gene_Dict[gene2])
        chr1 = bin_chrs[bin1]
        chr2 = bin_chrs[bin2]
        #skip if pair is on the same chromosome
        if chr1==chr2:
            if verbose:
                print "%s and %s are on the same chromosome. They are %r and %r" %(gene1,gene2,chr1,chr2)
            continue
        if bin1<bin2:
            score=trans_contacts.get((bin1,bin2),0)
        if bin2<bin1:
            score=trans_contacts.get((bin2,bin1),0)
        scoreList.append([gene1,gene2,score])
        colocalizationList.append(score)
        #if gene1=='TTN' or gene2=='TTN': #For troubleshooting of the score counts
        #    print "%s and %s are in bins %r and %r and have a socre of %r." %(gene1,gene2,bin1,bin2,score)
    #runs the permutations
    for i in range(perms):
        #creates a new random gene list for matched chromosomes
        randomGeneList = list()
        for gene in genes:
            gene_bin = int(gene_Dict[gene])
            gene_chr = bin_chrs[gene_bin]
            added=False
            #makes sure that the same gene wasn't already chosen - has to use this method so that it controls for the multiple chromosomes, can't just random sample all
            while added==False:
                newGene = random.sample(bkg_Dict[gene_chr],1)[0]
                if newGene not in randomGeneList:
                    randomGeneList.append(newGene)
                    added=True
        #iterates over the random gene list for each time
        for gene_pair in list(itertools.combinations(randomGeneList, 2)) :
            gene1 = gene_pair[0]
            gene2 = gene_pair[1]
            bin1 = int(gene_Dict[gene1])
            bin2 = int(gene_Dict[gene2])
            chr1 = bin_chrs[bin1]
            chr2 = bin_chrs[bin2]
            #skip if pair is on the same chromosome
            if chr1==chr2:
                if verbose and i==1:
                    print "%s and %s are on the same chromosome. They are %r and %r" %(gene1,gene2,chr1,chr2)
                continue
            if bin1<bin2:
                permutation_Dict[i].append(trans_contacts.get((bin1,bin2),0))
            if bin2<bin1:
                permutation_Dict[i].append(trans_contacts.get((bin2,bin1),0))
    return colocalizationList,scoreList,permutation_Dict

permutations = 1000
res=500000





#Load configuration File
c = ConfigParser.ConfigParser()
c.read(sys.argv[1])
sec = c.sections()[0]

#Assign Variables
gene_File = c.get(sec,'genes_file')
gene_List = c.get(sec,'gene_list')
contact_File = c.get(sec,'contact_file')
bed_File = c.get(sec,'bed_file')
perm_output_File = c.get(sec,'perm_output')
score_output_File = c.get(sec,'score_output')
bkg_File = c.get(sec,'background')


#opens the two other files to read and the new output file
all_Genes = aopen.open(gene_File,'r')
contacts = aopen.open(contact_File,'r')
bin_bed = aopen.open(bed_File,'r')
bkg_genes = aopen.open(bkg_File,'r')
output = aopen.open(perm_output_File,'w')
score_output = aopen.open(score_output_File,'w')

#creates a dictionary for the bin information based on the lower number of the bed file
#creates a dictionary of bins and their corresponding chromosomes
bin_data = {} #dictionary for chromosome location to get bin
bin_chrs = {} #dictionary for bin to get chromosome
i = 0
for line in bin_bed:
    bin_info = line.rstrip("\n").split("\t")
    bin_chr = bin_info[0]
    start = int(bin_info[1])
    bin_no = int(bin_info[3])
    bin_data[bin_chr,start]=bin_no
    bin_chrs[bin_no] = bin_chr
    i+=1

print "there are %i bins read in" %i

#creates a dictionary of chromosomes, within each dicitonary is a dictionary of the corresponding genes and their bin
chr_Dict = {}  #dictionary of chromosome with a set of genes
gene_Dict = {} #dictionary of genes and their bin
i=0
for line in all_Genes:
    gene_data = line.rstrip("\n").split("\t")
    gene_Name = str(gene_data[3])
    gene_Chr = str(gene_data[0])
    if gene_Chr.startswith('chr'): #tests if annotated with chr1 or 1 for chromosome name
        chromosome = gene_Chr[3:]
    else:
        chromosome = gene_Chr
    strand = str(gene_data[5]) #adjusts for which strand the chromosome is on to get start site
    if strand == '+':
        start_site = int(gene_data[1])
    else:
        start_site = int(gene_data[2])
    nearest_bin_start = (start_site / res)*res
    gene_bin = bin_data[chromosome,nearest_bin_start]
    #adds gene to dictionary of sets
    if chromosome in chr_Dict:
        chr_Dict[chromosome].add(gene_Name)
    else:
        chr_Dict[chromosome] = set()
        chr_Dict[chromosome].add(gene_Name)
    #adds gene and bin to their dictionary
    gene_Dict[gene_Name]=gene_bin
    i+=1
        
print "Gene information read in. There are %i number of genes" % i

#make a dictionary of the background gene list
bkg_Dict = {}
i=0
for line in bkg_genes:
    gene = line.rstrip("\n")
    gene_bin = gene_Dict.get(gene)
    chromosome = bin_chrs.get(gene_bin)
    if chromosome in bkg_Dict:
        bkg_Dict[chromosome].add(gene)
    else:
        bkg_Dict[chromosome] = set()
        bkg_Dict[chromosome].add(gene)
    i+=1

print "Background list read in. There are %i gene in the backgound set" %i

#genes to compare
genes = list()
for item in gene_List.split(','):
    if item in gene_Dict:
        genes.append(item)

#creates a dictionary of all trans contact information from ICED matrix
trans_contacts = {} # dictionary of trans contacts
i = 0
for line in contacts:
    contact_data = line.rstrip("\n").split("\t")
    bin1 = int(contact_data[0])    
    bin2 = int(contact_data[1])
    counts = float(contact_data[2])
    #check if a read is a cis or trans contact and adds counts to dictionary
    if bin_chrs[bin1] != bin_chrs[bin2]:
        trans_contacts[bin1,bin2] = counts
        i+=1

print "Contacts read in. There are %i trans contacts" %i


#calculates the co-localization score and does permutations tests
verbose=True
colocalizationList,scoreList,permutation_Dicts=getCoLocalizationData(genes,bkg_Dict,gene_Dict,bin_chrs,permutations,verbose)

for i in range(len(scoreList)):
    item=scoreList[i]
    newLine = str(item[0]) + "\t" + str(item[1]) + "\t" + str(item[2]) + "\n" 
    score_output.write(newLine)
score_output.close()

#write output from experimental data
experimentalData = sum(colocalizationList)
experimentalDataLine = "Experimental" + "\t" + str(experimentalData) + "\n"
print experimentalDataLine.rstrip("\n")
print len(colocalizationList)
output.write(experimentalDataLine)
   
j = 0
for i in range(permutations):
    randomTotal = sum(permutation_Dicts[i])
    if experimentalData>randomTotal: #figures out if our pair is greater than the random total
        j+=1
    randomDataLine = "Random" + "\t" + str(randomTotal) + "\n"
    output.write(randomDataLine)
    if j<10:
        print randomDataLine.rstrip("\n")
        print len(permutation_Dicts[i])

print "The experimental data ranks %i out of 1000 permutations" %j

output.close()

#Runs through the permutations test again for removing one less gene, mostly for curiosity to see if one gene is an outlier for the set
for oneLessGene in list(itertools.combinations(genes, len(genes)-1)):
    verbose=False
    colocalizationList,scoreList,permutation_Dicts=getCoLocalizationData(oneLessGene,bkg_Dict,gene_Dict,bin_chrs,permutations,verbose)

    experimentalData = sum(colocalizationList)
    experimentalDataLine = "Experimental" + "\t" + str(experimentalData) + "\n"
    j = 0
    for i in range(permutations):
        randomTotal = sum(permutation_Dicts[i])
        if experimentalData>randomTotal: #figures out if our pair is greater than the random total
            j+=1
        randomDataLine = "Random" + "\t" + str(randomTotal) + "\n"

    print "The experimental data ranks %i out of 1000 permutations, without gene %r" %(j,(set(genes)-set(oneLessGene)))


