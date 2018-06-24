#Python Script designed to calculate the relative contact probablity in cis for compartment information from HiC
#A-A, A-B, B-B interaction over all chromosomes
#Normalized to total contacts at a given distance
#Must have compartment information and interaction score at the same resolution
#Inputs given by a config file

#Written by Paul Fields
#Last Modified June 24, 2018

import os
import sys
import ConfigParser
import numpy as np

#Load configuration File
c = ConfigParser.ConfigParser()
c.read(sys.argv[1])
sec = c.sections()[0]

#Assign Variables
chr_File = c.get(sec,'chr_sizes') # Takes a file of chromosome and sizes (2 columns)
bin_File = c.get(sec,'bin_bed') # Takes a 4 col file (chr, start, end, bin)
contact_File = c.get(sec,'contact_file') # Takes a 3 col HiC-Pro output (bin1, bin2, score)
compartment_Calls = c.get(sec,'compartment_file') # Takes a BED4 format (chr, start, end, score/PC1)
output_File = c.get(sec,'output') # Creates an output file
resolution = int(c.get(sec,'resolution')) #Resolution from mapping

#opens the two other files to read and the new output file
chroms = open(chr_File,'r')
bins = open(bin_File,'r')
contacts = open(contact_File,'r')
compartments = open(compartment_Calls,'r')
output = open(output_File,'w')

#Creates a list of all the chromosomes of interest (only looking at knonw autosomes):
#Creates a dictionary of chromosomes and chromosome sizes
chrs=[]
chr_Sizes={}
for line in chroms:
    chrom=line.rstrip("\n").split("\t")
    if chrom[0].isdigit():
        chrs.append(chrom[0])
        chr_Sizes[chrom[0]]=int(chrom[1])

#Creates a dictionary both for bin to get coordinates and the reverse to give it coordinates and get the matching bin.
bin_data = {} #dictionary for chromosome location to get bin
bin_chrs = {} #dictionary for bin to get chromosome
for line in bins:
    bin_info = line.rstrip("\n").split("\t")
    bin_chr = bin_info[0]
    start = int(bin_info[1])
    bin_no = int(bin_info[3])
    bin_data[bin_chr,start]=bin_no
    bin_chrs[bin_no] = [bin_chr,start]

#Creates a dictionary of the bin with the compartment information
comp={}
for line in compartments:
    info = line.rstrip("\n").split("\t")
    if len(info)>1:
        bin_no = bin_data[info[0],int(info[1])]
        comp[bin_no]=float(info[3])

#Assuming Chr1 is Longest, sets up counting tables for reads
distCount_A = np.ones( ( chr_Sizes['1'] / resolution ) + 1 ); binCount_A = np.zeros( ( chr_Sizes['1'] / resolution ) + 1 )
distCount_B = np.ones( ( chr_Sizes['1'] / resolution ) + 1 ); binCount_B = np.zeros( ( chr_Sizes['1'] / resolution ) + 1 )
distCount_AB = np.ones( ( chr_Sizes['1'] / resolution ) + 1 ); binCount_AB = np.zeros( ( chr_Sizes['1'] / resolution ) + 1 )
distCount_Total = np.ones( ( chr_Sizes['1'] / resolution ) + 1 ); binCount_Total = np.zeros( ( chr_Sizes['1'] / resolution ) + 1 )

for line in contacts:
    contact_data=line.rstrip("\n").split("\t")
    bin1 = int(contact_data[0]) #Extracts the info from the triplet contact matrix    
    bin2 = int(contact_data[1])
    counts = float(contact_data[2])
    chr1 = bin_chrs[bin1][0] #Gets the chromosome for the two points
    chr2 = bin_chrs[bin2][0]
    comp1 = comp.get(bin1,0) #Gets the compartment value
    comp2 = comp.get(bin2,0)
    if chr1 == chr2: #Checks if they are on the same chromosome
        dist = bin2-bin1 #calculates the distance as a multiple of the resolution
        if ( comp1 == 0 or comp2 == 0 or bin1 == bin2):
            continue
        elif ( comp1 > 0 and comp2 > 0 ) :
            binCount_A[dist] += 1; distCount_A[dist] += counts
            binCount_Total[dist] += 1; distCount_Total[dist] += counts
        elif ( comp1 < 0 and comp2 < 0 ) :
            binCount_B[dist] += 1; distCount_B[dist] += counts
            binCount_Total[dist] += 1; distCount_Total[dist] += counts
        else :
            binCount_AB[dist] += 1; distCount_AB[dist] += counts
            binCount_Total[dist] += 1; distCount_Total[dist] += counts


 
binCount_A[np.where(binCount_A==0)] = 1; scaledCount_A = distCount_A/binCount_A
binCount_B[np.where(binCount_B==0)] = 1; scaledCount_B = distCount_B/binCount_B
binCount_AB[np.where(binCount_AB==0)] = 1; scaledCount_AB = distCount_AB/binCount_AB
binCount_Total[np.where(binCount_Total==0)] = 1; scaledCount_Total = distCount_Total/binCount_Total

relative_A=scaledCount_A/scaledCount_Total
relative_B=scaledCount_B/scaledCount_Total
relative_AB=scaledCount_AB/scaledCount_Total

header="Distance\tA_Compartment\tB_Compartment\tAB_Compartment\n"
output.write(header)
for i in range(len(distCount_A)):
     line="%i\t%f\t%f\t%f\n" % (i,relative_A[i],relative_B[i],relative_AB[i])
     output.write(line)

chroms.close()
bins.close()
contacts.close()
compartments.close()
output.close()
