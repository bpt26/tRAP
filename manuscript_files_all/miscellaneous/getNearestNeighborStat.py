#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import random
import numpy
import gzip
import math

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg

##########################
##### MAIN FUNCTIONS #####
##########################

def getNearestNeighborStat():
    for species in ['Btaur8','Chirc1','Hsapi38']:
        startTotRNA = {}
        endTotRNA = {}
        chromToStarts = {}
        chromToEnds = {}
        tRNAToActivity = {}
        mytRNAs = {}
        if not species == 'Hsapi38':
            for line in open(species+'tRNAClassificationsNewNoSegDups.txt'):
                splitLine = (line.strip()).split('\t')
                tRNAToActivity[splitLine[0]] = splitLine[2]
            for line in open(species+'tRNAHiConf.bed'):
                splitLine = (line.strip()).split('\t')
                if not splitLine[0] in chromToStarts:
                    chromToStarts[splitLine[0]] = []
                    chromToEnds[splitLine[0]] = []
                if splitLine[3] in tRNAToActivity:
                    chromToStarts[splitLine[0]].append(int(splitLine[1]))
                    chromToEnds[splitLine[0]].append(int(splitLine[2]))
                    startTotRNA[splitLine[0]+'_'+str(splitLine[1])] = splitLine[3]
                    endTotRNA[splitLine[0]+'_'+str(splitLine[2])] = splitLine[3]
        else:
            for line in open('humanCpGTrainingSet.tsv'):
                splitLine = (line.strip()).split('\t')
                tRNAToActivity[splitLine[0]] = splitLine[-1]
            for line in open(species+'tRNAHiConf.bed'):
                splitLine = (line.strip()).split('\t')
                if not splitLine[0] in chromToStarts:
                    chromToStarts[splitLine[0]] = []
                    chromToEnds[splitLine[0]] = []
                if splitLine[3] in tRNAToActivity:
                    chromToStarts[splitLine[0]].append(int(splitLine[1]))
                    chromToEnds[splitLine[0]].append(int(splitLine[2]))
                    startTotRNA[splitLine[0]+'_'+str(splitLine[1])] = splitLine[3]
                    endTotRNA[splitLine[0]+'_'+str(splitLine[2])] = splitLine[3]
        same = 0
        diff = 0
        for chrom in chromToStarts:
            myLen = len(chromToStarts[chrom])
            i = 0
            while i < myLen:
                mytRNA = startTotRNA[chrom+'_'+str((chromToStarts[chrom])[i])]
                myActivity = tRNAToActivity[mytRNA]
                if myLen > 1:
                    if i == 0:
                        nearest = startTotRNA[chrom+'_'+str((chromToStarts[chrom])[i+1])]
                    elif i == myLen-1:
                        nearest = startTotRNA[chrom+'_'+str((chromToStarts[chrom])[i-1])]
                    else:
                        test1 = (chromToStarts[chrom][i+1])-(chromToEnds[chrom][i])
                        test2 = (chromToStarts[chrom][i])-(chromToEnds[chrom][i-1])
                        if test1 < test2:
                            nearest = startTotRNA[chrom+'_'+str((chromToStarts[chrom])[i+1])]
                        else:
                            nearest = startTotRNA[chrom+'_'+str((chromToStarts[chrom])[i-1])]

                    if tRNAToActivity[nearest] == myActivity:
                        same += 1
                    else:
                        diff += 1
                i += 1
        print(same/float(same+diff)*100.0)
                






def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ' '.join(newList)

def main():
    getNearestNeighborStat()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit














