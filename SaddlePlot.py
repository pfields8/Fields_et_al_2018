#Distance Probablities By Distance - Generates Binned Decay curves - With Given Parameters
#Modified script from cooltools https://github.com/mirnylab/cooltools/tree/master/cooltools
#Can generate saddle plots for either cis or trans reads using triplet matrix with compartment scores
#Inputs are taken in as parameters

#Written by Paul Fields
#Last Modified June 24, 2018


#========================================================
#================MODULES=================================
#========================================================

import os
import sys
import string
import subprocess
import re
import numpy as np
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt


#========================================================
#================FUNCTIONS===============================
#========================================================

def chrSizes(chrFile):
    #Creates a list of all the chromosomes of interest (only looking at knonw autosomes):
    #Creates a dictionary of chromosomes and chromosome sizes
    chroms=open(chrFile,'r')
    chrs=[]
    chr_Sizes={}
    for line in chroms:
        chrom=line.rstrip("\n").split("\t")
        if chrom[0].isdigit():
            chrs.append(int(chrom[0]))
            chr_Sizes[int(chrom[0])]=int(chrom[1])
    chroms.close() 
    return chrs,chr_Sizes

def binInformation(binFile):
    #Creates a dictionary both for bin to get coordinates and the reverse to give it coordinates and get the matching bin.
    bins=open(binFile,'r')
    bin_data = {} #dictionary for chromosome location to get bin
    bin_chrs = {} #dictionary for bin to get chromosome
    for line in bins:
        bin_info = line.rstrip("\n").split("\t")
        bin_chr = bin_info[0]
        start = int(bin_info[1])
        bin_no = int(bin_info[3])
        bin_data[bin_chr,start]=bin_no
        bin_chrs[bin_no] = [bin_chr,start]
    bins.close()
    return bin_data,bin_chrs

def compartmentInfo(compartmentFile,bin_data):
    #Creates a dictionary of the bin with the compartment information
    compartments=open(compartmentFile,'r')
    comp={}
    for line in compartments:
        info = line.rstrip("\n").split("\t")
        if len(info)>1:
            bin_no = bin_data[info[0],int(info[1])]
            comp[bin_no]=float(info[3])
    compartments.close()
    return comp

def digitizeComps(compartments,segments):
    #Create an array for the compartment values given the specified number of bins 
    #Assigned by percentile
    n=segments
    temp_values=list(compartments.values())
    AB_values=sorted(temp_values, key=float)
    perc_edges = np.linspace(0, 100, n + 1)
    edges = np.percentile(AB_values, perc_edges)
    digitizedComparts = {}
    for key in compartments:
        digitizedComparts[key]=np.digitize(compartments[key],edges, right=True)
    return digitizedComparts,edges

def transSaddle(contactFile,bin_chrs,bin_data,digitizedComparts,segments):
    n_bins = segments

    interaction_sum   = np.zeros((n_bins, n_bins))
    interaction_count = np.ones((n_bins, n_bins))    

    contacts = open(contactFile,'r')
    for line in contacts:
        contact_data=line.rstrip("\n").split("\t")
        bin1 = int(contact_data[0]) #Extracts the info from the triplet contact matrix    
        bin2 = int(contact_data[1])
        counts = float(contact_data[2])
        chr1 = bin_chrs[bin1][0] #Gets the chromosome for the two points
        chr2 = bin_chrs[bin2][0]
        comp1 = digitizedComparts.get(bin1,0) #Gets the compartment value
        comp2 = digitizedComparts.get(bin2,0)
        if chr1 != chr2: #Checks if they are on the same chromosome
            if ( comp1 == 0 or comp2 == 0 ): #skips un-annotated compartments
                continue
            else:
                interaction_sum[comp1-1,comp2-1]   += 1
                interaction_count[comp1-1,comp2-1] += counts

    interaction_sum   += interaction_sum.T
    interaction_count += interaction_count.T

    total_sum = np.sum(interaction_sum)
    total_counts = np.sum(interaction_count)
    average_counts=float(total_counts/total_sum)
  
    #returns of the average interaction count for a given set of bins
    saddle_temp=interaction_count/interaction_sum
    saddle=saddle_temp/average_counts
    return saddle

def cisReads(contactFile,bin_chrs):
    #Resolution is set at 500000
    #Maximum distance is 100MB - length is 
    distCount = np.ones( ( 100000000 / 500000 ) + 1 ); binCount = np.zeros( (100000000 / 500000 ) + 1 )
    cis_reads={}

    contacts = open(contactFile,'r')
    for line in contacts:
        contact_data=line.rstrip("\n").split("\t")
        bin1 = int(contact_data[0]) #Extracts the info from the triplet contact matrix    
        bin2 = int(contact_data[1])
        counts = float(contact_data[2])
        chr1 = bin_chrs[bin1][0] #Gets the chromosome for the two points
        chr2 = bin_chrs[bin2][0]
        if chr1 == chr2: #Checks if they are on the same chromosome
            dist = bin2-bin1 #calculates the distance as a multiple of the resolution
            if dist>200: #Sets max distance for calculating average values at 100MB
                dist=200
            if ( bin1 == bin2):
                continue
            else:
                binCount[dist] += 1; distCount[dist] += counts
                cis_reads[bin1,bin2]=counts     
     
    #Can't dvidie by zero
    binCount[np.where(binCount==0)] = 1; averageCount = distCount/binCount
    return averageCount,cis_reads

def cisSaddle(contactFile,bin_chrs,bin_data,digitizedComparts,segments):
    n_bins = segments

    interaction_sum   = np.zeros((n_bins, n_bins))
    interaction_count = np.ones((n_bins, n_bins))

    averageCount,cis_reads=cisReads(contactFile,bin_chrs)

    for key in cis_reads:
        bin1=key[0]
        bin2=key[1]
        counts=cis_reads[key]
        comp1 = digitizedComparts.get(bin1,0) #Gets the compartment value
        comp2 = digitizedComparts.get(bin2,0)
        if ( comp1 == 0 or comp2 == 0 ): #skips un-annotated compartments
            continue
        else:
            dist = bin2-bin1 #calculates the distance as a multiple of the resolution
            if dist>200: #Sets max distance for calculating average values at 100MB
                dist=200
            average_value=averageCount[dist]
            normalized=counts/average_value
            interaction_sum[comp1-1,comp2-1]   += 1
            interaction_count[comp1-1,comp2-1] += normalized        

    interaction_sum   += interaction_sum.T
    interaction_count += interaction_count.T

    #returns of the average interaction count for a given set of bins
    saddle=interaction_count/interaction_sum
    return saddle


def make_saddle(binFile,contactFile,compartmentFile,contact_type,segments):
    #Create a matrix of average interactions given a digitized compartments
    #March 29, test run only designed for trans

    #Creates a dictionary both for bin to get coordinates and the reverse to give it coordinates and get the matching bin.
    bin_data,bin_chrs=binInformation(binFile)

    #Creates a dictionary of the bin with the compartment information
    comp=compartmentInfo(compartmentFile,bin_data)

    #Generates the digizted dictionary of compartments
    digitizedComparts,edges=digitizeComps(comp,segments)

    if contact_type=='trans':
        saddle=transSaddle(contactFile,bin_chrs,bin_data,digitizedComparts,segments)
    if contact_type=='cis':
        saddle=cisSaddle(contactFile,bin_chrs,bin_data,digitizedComparts,segments)

    return saddle,edges


def main():
    '''
    Main wrapper function
    '''
    from optparse import OptionParser

    usage = "This script takes triplet Hi-C contacts and generates a saddle plot for A B compartments, either in Trans or Cis. Usage: %prog [options] -r [Resolution]  -c [Contact File] -b [Bin File] -ab [Compartment Calls] -o [OUTPUTFILE]"

    parser = OptionParser(usage = usage)
    #required flags
    parser.add_option("-r","--resolution", dest="resolution",nargs = 1, default=500000,
                      help = "Set the resolution of the Hi-C Data")
    parser.add_option("-i","--contacts", dest="contacts",nargs = 1, default=None,
                      help = "Enter the triplet contact file")
    parser.add_option("-b","--bins", dest="bins",nargs = 1, default=None,
                      help = "Enter a 4 column bed file corresponding to cordinates and bin number.")
    parser.add_option("-s","--chr_sizes", dest="chr_sizes",nargs = 1, default=None,
                      help = "Enter a file with chromosome sizes for those chromosomes to be processed.")

    #output flag
    parser.add_option("-o","--output", dest="output",nargs = 1, default=None,
                      help = "Enter a name for the output file.")

    #additional options
    parser.add_option("-c","--compartments", dest="compartments",nargs = 1, default=None,
                      help = "Enter a bed file with a 4th column of the PC1 scores for compartment calls")
    parser.add_option("-t","--contact_type", dest="contact_type",nargs = 1, default=None,
                      help = "Generate a saddle plot for cis or trans data.")
    parser.add_option("-d","--segments", dest="segments",nargs = 1, default=None,
                      help = "Number of segments to generate for a matrix.")


    (options,args) = parser.parse_args()
    print('Options and arguments passed to script:')
    print(options)
    print(args)

    #Checks for resolution
    try:
        int(options.segments)
    except:
        print('Warning: No segments specified under -d flag. Using default of 10.')
        options.segments = 10
    segments = int(options.segments)

    #Checks for number of segments to make for bins
    try:
        int(options.resolution)
    except:
        print('Warning: No resolution specified under -r flag. Using default of 500KB.')
        options.resolution = 500000
    resolution = int(options.resolution)

    contact_type = options.contact_type
    if contact_type not in ['cis', 'trans']:
        raise ValueError("The allowed values for the contact_type "
                         "argument are 'cis' or 'trans'.")
    
    if options.contacts and options.bins and options.chr_sizes and options.compartments:
        if options.output == None:
            output = os.getcwd() + contactsFile.split('/')[-1]+'_saddle.txt'
            print('No output specified. Writing to %s' % (output))
        else:
            output = options.output

        contactFile = options.contacts
        binFile = options.bins
        compartmentFile = options.compartments

        saddle,edges=make_saddle(binFile,contactFile,compartmentFile,contact_type,segments)
        
        headerText=""
        for n in range(edges.size-2):
            headerText=headerText+"%f-%f\t" %(edges[n],edges[n+1])
        headerText=headerText+"%f-%f" %(edges[n],edges[n+1])

        np.savetxt(output,saddle,delimiter="\t",header=headerText)

    else:
        parser.print_help()

if __name__ == "__main__":
    main()


