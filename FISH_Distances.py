#Python Script for integrating the output from IMARIS for use in downstream applications
#Takes in the data for FISH probe locations, and distance to each other, and also distance to nuclei border
#Takes in nuclei area, vol, and coordinates of center of nuclei shape
#Merges multiple sets of images to a single output with cells that contain two loci each for a red and green probe.
#Inputs given in comma seperated list, files are comma seperated output from IMARIS

#Written by Paul Fields
#Last modified June 24, 2018

import sys
import ConfigParser
import csv
from itertools import izip
import math

def getNuclei(x,y,z,nuclei_dict):
    dist = {}
    for key, value in nuclei_dict.iteritems():
         distance = math.sqrt( (x-value[0])**2 + (y-value[1])**2 + (z-value[2])**2 )
         dist[key] = distance
    min_key = min(dist, key=dist.get)
    return min_key

#New way to read in configuration file
c = ConfigParser.SafeConfigParser()
c.read(sys.argv[1])
sec = c.sections()[0]

#Assign Variables - set default variable
nucleiArea = c.get(sec,'nucleiArea').split(',')
nucleiVol = c.get(sec,'nucleiVol').split(',')
nucleiCoordinate = c.get(sec,'nucleiCoordinate').split(',')
GFP_coordinates = c.get(sec,'GFP_coordinates').split(',')
GFP_nucleiDist = c.get(sec,'GFP_nuclei').split(',')
GFP_dists = c.get(sec,'GFP_dists').split(',')
RFP_coordinates = c.get(sec,'RFP_coordinates').split(',')
RFP_nucleiDist = c.get(sec,'RFP_nuclei').split(',')
RFP_dists = c.get(sec,'RFP_dists').split(',')
outfile = c.get(sec,'Outfile')

#Get number of files
numFiles = len(nucleiArea)

#Creates dictionary of dictionary for nuclei data
nuclei_stats = {}
nuclei_GFP = {}
nuclei_RFP = {}
GFP_dict = {}
RFP_dict = {}
GFP_edgeDist = {}
RFP_edgeDist = {}
nucleiToUse = {}

for i in range(numFiles):
    nuclei_stats[i] = {}
    nuclei_GFP[i] = {}
    nuclei_RFP[i] = {}
    GFP_dict[i] = {}
    RFP_dict[i] = {}
    GFP_edgeDist[i] = {}
    RFP_edgeDist[i] = {}
    nucleiToUse[i] = []

#Read in nuclei area for all files
for i in range(numFiles):
    with open(nucleiArea[i],'r') as f1, open(nucleiVol[i],'r') as f2, open(nucleiCoordinate[i],'r') as f3:
        #reader 1 = Area, reader 2 = volume, reader 3 = coordinates
        reader1 = csv.reader(f1)
        reader2 = csv.reader(f2)
        reader3 = csv.reader(f3)
        for v1, v2, v3 in izip(reader1, reader2, reader3):
            if len(v1)==0:
                continue
            try: #Checks that its not a header row, that it can convert to floats
                ID=v1[4]
                area=float(v1[0])
                vol=float(v2[0])
                x=float(v3[0])
                y=float(v3[1])
                z=float(v3[2])
                nuclei_stats[i][ID]=[x,y,z,area,vol] #stores data as x,y,z,area,vol
            except:
                continue

#Analyze the GFP spots for all files
for i in range(numFiles):
    with open(GFP_coordinates[i],'r') as f4, open(GFP_dists[i],'r') as f5, open(GFP_nucleiDist[i],'r') as f8:
        r4 = csv.reader(f4)
        r5 = csv.reader(f5)
        r8 = csv.reader(f8)
        for v1,v2,v3 in izip(r4,r5,r8):
            if len(v1)==0:
                continue
            try: #Checks for headers and if can't fit coordinates skips row
                x1=float(v1[0])
                y1=float(v1[1])
                z1=float(v1[2])
                toRFP=float(v2[0]) #Min distance to an RFP spot
                toEdge=float(v3[0]) #Min distance to nuclei edge
            except:
                continue
            ID=v1[7]
            nuclei=getNuclei(x1,y1,z1,nuclei_stats[i])
            nuclei_GFP[i].setdefault(nuclei,[]).append(ID)
            GFP_dict[i][ID]=toRFP
            GFP_edgeDist[i][ID]=toEdge

#Analyze the RFP spots for all files
for i in range(numFiles):
    with open(RFP_coordinates[i],'r') as f6, open(RFP_dists[i],'r') as f7, open(RFP_nucleiDist[i],'r') as f9:
        r6 = csv.reader(f6)
        r7 = csv.reader(f7)
        r9 = csv.reader(f9)
        for v1,v2,v3 in izip(r6,r7,r9):
            if len(v1)==0:
                continue
            try: #Checks for headers and if can't fit coordinates skips row
                x1=float(v1[0])
                y1=float(v1[1])
                z1=float(v1[2])
                toGFP=float(v2[0]) #Min distance to a GFP spot
                toEdge=float(v3[0]) #Min distance to nuclei edge
            except:
                continue
            ID=v1[7]
            nuclei=getNuclei(x1,y1,z1,nuclei_stats[i])
            nuclei_RFP[i].setdefault(nuclei,[]).append(ID)
            RFP_dict[i][ID]=toGFP
            RFP_edgeDist[i][ID]=toEdge

#Generate list of spots to analyze for all files
for i in range(numFiles):
    for key in nuclei_stats[i].iterkeys():
        if len(nuclei_GFP[i].get(key,[]))==2 and len(nuclei_RFP[i].get(key,[]))==2:
            nucleiToUse[i].append(key)

#Output data
toWrite=open(outfile,'w')
header = "Nuclei_ID\tArea\tVol\tToRFP_1\tToRFP_2\tToGFP_1\tToGFP_2\tGFPtoEdge_1\tGFPtoEdge_2\tRFPtoEdge_1\tRFPtoEdge_2\n"
toWrite.write(header)

for i in range(numFiles):
    for nuc in nucleiToUse[i]: #Iterates through nuclei to use in each file
        area = nuclei_stats[i][nuc][3] #Gets area and volume
        vol = nuclei_stats[i][nuc][4]
        toRFP = [] #Distances for GFP and RFP to other type
        toGFP = []
        GFPtoEdge = [] #Distances for GFP and RFP to nuclei edge
        RFPtoEdge = []
        for j in range(2): #Each cell should have only two spots
            spotGFP = nuclei_GFP[i][nuc][j] #each are a list of two nuclei IDs corresponding to that nucleus
            spotRFP = nuclei_RFP[i][nuc][j]
            toRFP.append(GFP_dict[i][spotGFP]) #Gets the distances to the other spots
            toGFP.append(RFP_dict[i][spotRFP])
            GFPtoEdge.append(GFP_edgeDist[i][spotGFP])
            RFPtoEdge.append(RFP_edgeDist[i][spotRFP])        

        outline = "%r\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (int(nuc),area,vol,toRFP[0],toRFP[1],toGFP[0],toGFP[1],GFPtoEdge[0],GFPtoEdge[1],RFPtoEdge[0],RFPtoEdge[1]) #add new line for output writing
        toWrite.write(outline)
        print "Output from Files #%i" % (i+1)
        print outline 
        
toWrite.close()


