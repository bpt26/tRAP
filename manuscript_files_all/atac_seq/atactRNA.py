#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# redoGetFeatures.py

import sys
import os
import time
import random
import numpy
import gzip
import math

##########################
##### MAIN FUNCTIONS #####
##########################

"""
This script is used to map ATAC-seq data from Foissac et al 2018 to regions surrounding tRNA genes by 250 nt 
on either side. These files have peak data from several replicates. This script takes the mean for each tRNA, 
which can later be plotted versus our predictions.
"""

def atactRNA(abb, bed):

    chromConverter = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/'+abb+'/'+abb+'_assembly_report.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 8:
            if abb in ['Btaur8', 'Chirc1']:
                chromConverter[splitLine[4]] = str(splitLine[6])
                chromConverter[str(splitLine[2])] = str(splitLine[6])
            elif abb == 'Sscro11':
                chromConverter[str(splitLine[2])] = str(splitLine[4])
    sys.stderr.write("Finished reading assembly report.\n")

    coordTotRNA = {}
    myChroms = {}

    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/'+abb+'/'+abb+'_halChroms_hiConftRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        myChrom = str(splitLine[0])
        if not myChrom in myChroms:
            myChroms[myChrom] = True
        for k in range(int(splitLine[1])-250,int(splitLine[2])+250):
            coordTotRNA[myChrom+'-'+str(k)] = splitLine[3]
    sys.stderr.write("Finished reading tRNA bed.\n")

    nameMap = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/'+abb+'/'+abb+'NameMap.txt'):
        splitLine = (line.strip()).split('\t')
        nameMap[splitLine[0]] = splitLine[1]
        nameMap[splitLine[1]] = splitLine[0]


    tRNAToPeaks = {}
    myOut = '/public/groups/corbettlab/tRNA/classifier/ATAC_SEQ/BEDS/'+abb+'_ATAC_SEQ.out'
    for line in gzip.open('/public/groups/corbettlab/tRNA/classifier/ATAC_SEQ/BEDS/'+bed):
        splitLine = (line.strip()).split('\t')
        if str(splitLine[0]) in myChroms or str(splitLine[0]) in chromConverter:
            if not str(splitLine[0]) in myChroms and str(splitLine[0]) in chromConverter:
                myChrom = chromConverter[str(splitLine[0])]
            mytRNA = ''
            for k in range(int(splitLine[1]),int(splitLine[2])):
                if myChrom+'-'+str(k) in coordTotRNA:
                    mytRNA = coordTotRNA[myChrom+'-'+str(k)]
            if not mytRNA == '':
                tRNAToPeaks[mytRNA] = takeMean(splitLine[6:])
    myString = ''
    for k in tRNAToPeaks:
        myString += abb+'\t'+k+'\t'+nameMap[k]+'\t'+str(tRNAToPeaks[k])+'\n'
    open(myOut, 'w').write(myString)
                

def takeMean(myList):
    avg = 0.0
    for k in myList:
        avg += float(k)
    return(avg/float(len(myList)))

def main():
    atactRNA('Sscro11', 'sus_scrofa_mergedpeaks.peaknb.allexp.readnb.bed.gz')
    atactRNA('Chirc1', 'capra_hircus_mergedpeaks.peaknb.allexp.readnb.bed.gz')
    atactRNA('Btaur8', 'bos_taurus_mergedpeaks.peaknb.allexp.readnb.bed.gz')
    

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

