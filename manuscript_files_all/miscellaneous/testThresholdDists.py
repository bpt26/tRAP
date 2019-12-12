#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 3/9/2018
# fixIsotypeSpecific.py

import sys
import os
import time
import random
import numpy
import gzip
import math
import scipy
from scipy.stats import mannwhitneyu

##########################
##### MAIN FUNCTIONS #####
##########################


def testCloseness():
    THRESHOLD_DIST = int(sys.argv[1])
    myDist = 0
    myActivity = 'null'
    mytRNA = 'null'
    myDistEqual = 0.0
    myDistAndActEqual = 0.0
    matchDists = []
    unmatchDists = []
    singleActivity = []
    for line in open('/Users/Bryan/Desktop/classifier/mm10-tRNAs/mouseCoordsTest.txt', 'U'):
        splitLine = (line.strip()).split('\t')
        if splitLine[-1] != 'x':
            myPrevDist = myDist
            myPrevActivity = myActivity
            myPrevtRNA = mytRNA
            myDist = int(splitLine[-1])
            myActivity = splitLine[2]
            mytRNA = splitLine[0]
            if myDist == myPrevDist:
                myDistEqual += 1.0
                if myActivity == myPrevActivity:
                    myDistAndActEqual += 1.0
                    matchDists.append(myDist)
                else:
                    unmatchDists.append(myDist)
            if myDist >= THRESHOLD_DIST:
                singleActivity.append(myActivity)

    ##### WRITE SOMETHING THAT WILL FIND THE NEAREST PROTEIN CODING GENE 
    # FOR EACH TRNA AND LOOK AT THE DISTRIBUTION OF ACTIVITY GROUPS FOR THIS

    print(myDistEqual)
    print(myDistAndActEqual)
    print(sorted(matchDists))
    print(sorted(unmatchDists))
    print(numpy.mean(matchDists))
    print(numpy.mean(unmatchDists))
    print(numpy.median(matchDists))
    print(numpy.median(unmatchDists))
    #print(singleActivity.count('active'))
    #print(singleActivity.count('active')/float(len(singleActivity)))
    print(singleActivity.count('A'))
    print(singleActivity.count('A')/float(len(singleActivity)))
    print(singleActivity.count('B'))
    print(singleActivity.count('B')/float(len(singleActivity)))
    print(singleActivity.count('C'))
    print(singleActivity.count('C')/float(len(singleActivity)))
    print(singleActivity.count('D'))
    print(singleActivity.count('D')/float(len(singleActivity)))
    print(singleActivity.count('E'))
    print(singleActivity.count('E')/float(len(singleActivity)))


def gettRNAsWithinDist():

    THRESHOLD_DIST = int(sys.argv[1])
    tRNAStartDict = {}
    tRNAEndDict = {}
    startTotRNA = {}
    tRNAToStart = {}
    endTotRNA = {}
    tRNAToEnd = {}
    tRNAToActivity = {}

    for line in open('/Users/Bryan/Desktop/classifier/mouse_sup_table.csv'):
        splitLine = (line.strip()).split(',')
        if not splitLine[0] == 'tRNA':
            mytRNA = splitLine[0]
            if splitLine[1] in ['A','B','C']:
                tRNAToActivity[mytRNA] = 'active'
            else:
                tRNAToActivity[mytRNA] = 'inactive'


    for k in range(1,23):
        tRNAStartDict['chr'+str(k)] = []
        tRNAEndDict['chr'+str(k)] = []
    tRNAStartDict['chrX'] = []
    tRNAEndDict['chrX'] = []
    for line in open('/Users/Bryan/Desktop/classifier/mm10-tRNAs/mm10-tRNAs.bed', 'U'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            mytRNA = splitLine[3]
            myChrom = str(splitLine[0])
            (tRNAStartDict[myChrom]).append(int(splitLine[1]))
            (tRNAEndDict[myChrom]).append(int(splitLine[2]))
            tRNAToStart[mytRNA] = myChrom+'-'+str(splitLine[1])
            tRNAToEnd[mytRNA] = myChrom+'-'+str(splitLine[2])
            startTotRNA[myChrom+'-'+str(splitLine[1])] = mytRNA
            endTotRNA[myChrom+'-'+str(splitLine[2])] = mytRNA
    for chrom in tRNAStartDict:
        tRNAStartDict[chrom] = sorted(tRNAStartDict[chrom])
        tRNAEndDict[chrom] = sorted(tRNAEndDict[chrom])

    tRNAsInRange = {}
    for chrom in tRNAStartDict:
        for start in tRNAStartDict[chrom]:
            mytRNA = startTotRNA[chrom+'-'+str(start)]
            tRNAsInRange[mytRNA] = []
            for otherStart in tRNAStartDict[chrom]:
                if otherStart >= start-THRESHOLD_DIST and otherStart <= start+THRESHOLD_DIST+1:
                    myOthertRNA = startTotRNA[chrom+'-'+str(otherStart)]
                    if (myOthertRNA not in tRNAsInRange[mytRNA]) and (myOthertRNA != mytRNA):
                        (tRNAsInRange[mytRNA]).append(myOthertRNA)
            for end in tRNAEndDict[chrom]:
                if end >= start-THRESHOLD_DIST and end <= start+THRESHOLD_DIST+1:
                    myOthertRNA = endTotRNA[chrom+'-'+str(end)]
                    if (myOthertRNA not in tRNAsInRange[mytRNA]) and (myOthertRNA != mytRNA):
                        (tRNAsInRange[mytRNA]).append(myOthertRNA)

    tRNAList = []
    myOut = '/Users/Bryan/Desktop/classifier/MousetRNAsWithin'+str(THRESHOLD_DIST)+'.txt'
    open(myOut, 'w').write('tRNA\tcount\n')
    for tRNA in tRNAsInRange:
        tRNAList.append(tRNA)
    tRNAList = sorted(tRNAList, key=lambda k: len(tRNAsInRange[k]))
    tRNAZeroActivity = []
    for tRNA in tRNAList:
        if len(tRNAsInRange[tRNA]) == 0 and tRNA in tRNAToActivity:
            tRNAZeroActivity.append(tRNAToActivity[tRNA])
        #sys.stderr.write(tRNA+': '+str(len(tRNAsInRange[tRNA]))+'\t'+tRNAToActivity[tRNA]+'\n')
        if tRNA in tRNAToActivity:
            open(myOut, 'a').write(tRNA+'\t'+str(len(tRNAsInRange[tRNA]))+'\t'+str(tRNAToActivity[tRNA])+'\n')
    sys.stderr.write('active: '+str(tRNAZeroActivity.count('active'))+'\t')
    sys.stderr.write('inactive: '+str(tRNAZeroActivity.count('inactive'))+'\n')

    myActivetRNAs = []
    myInactivetRNAs = []
    for line in open(myOut):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            if splitLine[2] == 'active':
                myActivetRNAs.append(int(splitLine[1]))
            else:
                myInactivetRNAs.append(int(splitLine[1]))

    print(str(THRESHOLD_DIST)+'\ttRNA')
    print(mannwhitneyu(myActivetRNAs,myInactivetRNAs))




def getProteinCodingGenesWithinDist():

    THRESHOLD_DIST = int(sys.argv[1])
    tRNAStartDict = {}
    tRNAEndDict = {}
    startTotRNA = {}
    tRNAToStart = {}
    endTotRNA = {}
    tRNAToEnd = {}
    tRNAToActivity = {}

    protStartDict = {}
    protEndDict = {}
    startToProt = {}
    protToStart = {}
    endToProt = {}
    protToEnd = {}

    for line in open('/Users/Bryan/Desktop/classifier/mouse_sup_table.csv'):
        splitLine = (line.strip()).split(',')
        if not splitLine[0] == 'tRNA':
            mytRNA = splitLine[0]
            if splitLine[1] in ['A','B','C']:
                tRNAToActivity[mytRNA] = 'active'
            else:
                tRNAToActivity[mytRNA] = 'inactive'


    for k in range(1,23):
        tRNAStartDict['chr'+str(k)] = []
        tRNAEndDict['chr'+str(k)] = []
    tRNAStartDict['chrX'] = []
    tRNAEndDict['chrX'] = []
    for line in open('/Users/Bryan/Desktop/classifier/mm10-tRNAs/mm10-tRNAs.bed', 'U'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            mytRNA = splitLine[3]
            myChrom = str(splitLine[0])
            (tRNAStartDict[myChrom]).append(int(splitLine[1]))
            (tRNAEndDict[myChrom]).append(int(splitLine[2]))
            tRNAToStart[mytRNA] = myChrom+'-'+str(splitLine[1])
            tRNAToEnd[mytRNA] = myChrom+'-'+str(splitLine[2])
            startTotRNA[myChrom+'-'+str(splitLine[1])] = mytRNA
            endTotRNA[myChrom+'-'+str(splitLine[2])] = mytRNA
    for chrom in tRNAStartDict:
        tRNAStartDict[chrom] = sorted(tRNAStartDict[chrom])
        tRNAEndDict[chrom] = sorted(tRNAEndDict[chrom])


    for k in range(1,23):
        protStartDict['chr'+str(k)] = []
        protEndDict['chr'+str(k)] = []
    protStartDict['chrX'] = []
    protEndDict['chrX'] = []
    for line in open('/Users/Bryan/Desktop/classifier/mouseKnownGene.txt', 'U'):
        splitLine = (line.strip()).split('\t')
        myGene = splitLine[0]
        myChrom = str((splitLine[1]))
        if myChrom in protStartDict:
            (protStartDict[myChrom]).append(int(splitLine[3]))
            (protEndDict[myChrom]).append(int(splitLine[4]))
            protToStart[myGene] = myChrom+'-'+str(splitLine[3])
            protToEnd[myGene] = myChrom+'-'+str(splitLine[4])
            startToProt[myChrom+'-'+str(splitLine[3])] = myGene
            endToProt[myChrom+'-'+str(splitLine[4])] = myGene
    for chrom in protStartDict:
        protStartDict[chrom] = sorted(protStartDict[chrom])
        protEndDict[chrom] = sorted(protEndDict[chrom])


    protsInRange = {}
    for chrom in tRNAStartDict:
        for start in tRNAStartDict[chrom]:
            mytRNA = startTotRNA[chrom+'-'+str(start)]
            protsInRange[mytRNA] = []
            for otherStart in protStartDict[chrom]:
                if otherStart >= start-THRESHOLD_DIST and otherStart <= start+THRESHOLD_DIST+1:
                    myProt = startToProt[chrom+'-'+str(otherStart)]
                    if (myProt not in protsInRange[mytRNA]):
                        (protsInRange[mytRNA]).append(myProt)
            for end in protEndDict[chrom]:
                if end >= start-THRESHOLD_DIST and end <= start+THRESHOLD_DIST+1:
                    myProt = endToProt[chrom+'-'+str(end)]
                    if (myProt not in protsInRange[mytRNA]):
                        (protsInRange[mytRNA]).append(myProt)

    tRNAList = []
    myOut = '/Users/Bryan/Desktop/classifier/mouseProtsWithin'+str(THRESHOLD_DIST)+'.txt'
    open(myOut, 'w').write('tRNA\tcount\n')
    for tRNA in protsInRange:
        if tRNA in tRNAToActivity:
            tRNAList.append(tRNA)
    tRNAList = sorted(tRNAList, key=lambda k: (len(protsInRange[k]), tRNAToActivity[k]))
    tRNAZeroActivity = []
    for tRNA in tRNAList:
        if len(protsInRange[tRNA]) == 0:
            tRNAZeroActivity.append(tRNAToActivity[tRNA])
        #sys.stderr.write(tRNA+': '+str(len(protsInRange[tRNA]))+'\t'+tRNAToActivity[tRNA]+'\n')
        open(myOut, 'a').write(tRNA+'\t'+str(len(protsInRange[tRNA]))+'\t'+str(tRNAToActivity[tRNA])+'\n')
    sys.stderr.write('active: '+str(tRNAZeroActivity.count('active'))+'\t')
    sys.stderr.write('inactive: '+str(tRNAZeroActivity.count('inactive'))+'\n')

    myActivetRNAs = []
    myInactivetRNAs = []
    for line in open(myOut):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            if splitLine[2] == 'active':
                myActivetRNAs.append(int(splitLine[1]))
            else:
                myInactivetRNAs.append(int(splitLine[1]))

    print(str(THRESHOLD_DIST)+'\tprotein')
    print(mannwhitneyu(myActivetRNAs,myInactivetRNAs))



#testCloseness()
# gettRNAsWithinDist()
getProteinCodingGenesWithinDist()

