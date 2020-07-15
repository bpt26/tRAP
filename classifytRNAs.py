#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2019
# classifytRNAs.py

import sys
import os
import time
import random
import numpy
import gzip
import math
import argparse
import sklearn
import numpy as np
from sklearn.impute import SimpleImputer
from sklearn.datasets import make_classification
from sklearn.ensemble import RandomForestClassifier

"""
This program classifies tRNAs using a random forest classifier.
"""


class CommandLine(object):
    """Handles the input arguments from the command line. Manages 
    the argument parser.

    Methods:
    Other than initialization, no methods are present, as its purpose is 
    simply to handle what is passed into the command line and pass that 
    into the class that performs the searching algorithm."""

    def __init__(self, inOpts=None):
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line input using argparse.
        '''
        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("-b", "--inputBed", help="Input .bed"+
            " containing genomic coordinates of all high-confidence tRNA genes.")
        self.parser.add_argument("-e", "--inputBedExt", help="Input .bed containing genomic coordinates of"+
            " all high-confidence tRNA genes and 350 bases upstream and downstream where possible.")
        self.parser.add_argument("-c", "--phastConsElements", help="PhastCons elements"+
            " found in or near high-confidence tRNA genes.")
        self.parser.add_argument("-w", "--inputWig", help="Input .wig file containing PhyloP"+
            " scores for all bases in or up to 20 bases away from high-confidence tRNA genes.")
        self.parser.add_argument("-t", "--inputOut", help="Input .out file from tRNAscan-SE 2.0"+
            " containing bit scores of all high-confidence tRNA genes.")
        self.parser.add_argument("-m", "--inputMFE", help="Input .mfe"+
            " file containing minimum free-energy scores for all high-confidence tRNA genes.")
        self.parser.add_argument("-f", "--inputFa", help="Input .fa file containing sequences for the regions"+
            "spanning 350 bases upstream to 350 bases downstream of all high-confidence tRNA genes.")
        self.parser.add_argument("-g", "--inputGFF", help="Input gff file in .bed format containing all exons"+
            " annotated in the genome of interest.")
        self.parser.add_argument("-l", "--chromLengths", help="Optional but recommended: File containing lengths"+
            " for all chromosomes in genome of interest (each line should be chrom<tab>0<tab>length<end>.", default='')
        self.parser.add_argument("-d", "--trainingFile", help="File containing human data for training the classifier.", default='')
        self.parser.add_argument("-a", "--activityFile", help="File containing activity data for each tRNA gene.", default='')
        self.parser.add_argument("-s", "--segDups", help="List of segmental duplications to exclude from prediction.", default='')
        self.parser.add_argument("-o", "--outputFile", help="The path to"+
            " and the filename of the output predictions file you are creating.", default='')
        self.parser.add_argument("-x", "--simplified", help="Flag for using the simplified (no annotation- or "+
            "alignment-based features.", default=False)
        if inOpts is None:
            self.args = vars(self.parser.parse_args())
        else:
            self.args = vars(self.parser.parse_args(inOpts))
        self.args = self.parser.parse_args()

class fileConverter(object):

    def __init__(self, inputBed, inputBedExt, phastConsElements, inputWig, inputOut, inputMFE, inputFa, inputGFF, chromLengths, trainingFile, activityFile, segDups, outputFile, simplified):
        self.inputBed = inputBed
        self.inputBedExt = inputBedExt
        self.phastConsElements = phastConsElements
        self.inputWig = inputWig
        self.inputOut = inputOut
        self.inputMFE = inputMFE
        self.inputFa = inputFa
        self.inputGFF = inputGFF
        self.chromLengths = chromLengths
        self.trainingFile = trainingFile
        self.activityFile = activityFile
        self.segDups = segDups
        self.outputFile = outputFile
        self.simplified = simplified

    def getFeatures(self):

        ########################
        ########################
        ### LOAD tRNA COORDS ###
        ########################
        ########################

        coordTotRNA = {}
        coordToSprinzl = {}
        tRNAToCoords = {}
        myHaltRNAs = []
        ACtoCounts = {}
        tRNAToStart = {}
        tRNAToEnd = {}
        tRNAToStrand = {}
        chromTotRNAs = {}

        chromToLen = {}
        for line in open(self.chromLengths):
            splitLine = (line.strip()).split('\t')
            if len(splitLine) == 3:
                chromToLen[splitLine[0]] = int(splitLine[2])

        coordToCons = {}
        if self.phastConsElements != '':
            for line in open(self.phastConsElements):
                splitLine = (line.strip()).split('\t')
                myChrom = splitLine[0]
                if len(chromToLen) != 0 and not myChrom in chromToLen:
                    for i in range(1,10):
                        if myChrom+'.'+str(i) in chromTotRNAs: 
                            myChrom = myChrom+'.'+str(i)
                myStart = int(splitLine[1])
                myEnd = int(splitLine[2])
                myScore = int(splitLine[4])
                if myEnd-myStart <= 100 and myScore >= 200:
                    for k in range(myStart, myEnd+1):
                        coordToCons[myChrom+'__'+str(k)] = True

        for line in open(self.inputBed):
            splitLine = (line.strip()).split('\t')
            myChrom = splitLine[0]
            myStart = int(splitLine[1])
            myEnd = int(splitLine[2])
            mytRNA = splitLine[3]
            myStrand = splitLine[5]

            if not myChrom in chromTotRNAs:
                chromTotRNAs[myChrom] = []
            chromTotRNAs[myChrom].append(mytRNA)

            tRNAToStart[mytRNA] = myChrom+'__'+str(myStart)
            tRNAToEnd[mytRNA] = myChrom+'__'+str(myEnd)
            tRNAToStrand[mytRNA] = str(myStrand)

            myAC = mytRNA[-3:]
            if not myAC in ACtoCounts:
                ACtoCounts[myAC] = 0
            ACtoCounts[myAC] += 1

            tRNAToCoords[mytRNA] = []
            myHaltRNAs.append(mytRNA)
            if int(splitLine[9]) > 1:
                mySizes = (splitLine[10]).split(',')[:-1]
                myStarts = (splitLine[11]).split(',')[:-1]
                for i in range(0, len(mySizes)):
                    for k in range(myStart+int(myStarts[i]), myStart+int(myStarts[i])+int(mySizes[i])):
                        coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                        coordToSprinzl[myChrom+'__'+str(k)] = 1
                        (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))
            else:
                for k in range(myStart, myEnd):
                    coordToSprinzl[myChrom+'__'+str(k)] = 1
                    coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))

            if myChrom+'__'+str(myStart) in coordToCons:
                while (myChrom+'__'+str(myStart) in coordToCons and myStart-20 > 0):
                    coordTotRNA[myChrom+'__'+str(myStart)] = mytRNA
                    coordToSprinzl[myChrom+'__'+str(myStart)] = 1
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(myStart))
                    myStart -= 1

            if myChrom+'__'+str(myEnd) in coordToCons:
                if myChrom in chromToLen:
                    while (myChrom+'__'+str(myEnd) in coordToCons and myEnd+20 < chromToLen[myChrom]):
                        coordTotRNA[myChrom+'__'+str(myEnd)] = mytRNA
                        coordToSprinzl[myChrom+'__'+str(myEnd)] = 1
                        (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(myEnd))
                        myEnd += 1

            if myStrand == '+':
                for k in range(myStart-20, myStart):
                    coordToSprinzl[myChrom+'__'+str(k)] = 5
                    coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))
                for k in range(myEnd, myEnd+10):
                    coordToSprinzl[myChrom+'__'+str(k)] = 3
                    coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))

            elif myStrand == '-':
                for k in range(myStart-10, myStart):
                    coordToSprinzl[myChrom+'__'+str(k)] = 3
                    coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))
                for k in range(myEnd, myEnd+20):
                    coordToSprinzl[myChrom+'__'+str(k)] = 5
                    coordTotRNA[myChrom+'__'+str(k)] = mytRNA
                    (tRNAToCoords[mytRNA]).append(myChrom+'__'+str(k))

        #########################
        #########################
        ### GET PHYLOP SCORES ###
        #########################
        #########################

        if self.simplified == False:
            sys.stderr.write("Getting PhyloP features...\n")

            coordToPhyloP = {}
            for line in open(self.inputWig):
                splitLine = (line.strip()).split()
                if (line.strip()).startswith('fixedStep'):
                    myChrom = (splitLine[1]).split('=')[-1]
                    if not myChrom in chromTotRNAs:
                        for i in range(1,10):
                            if myChrom+'.'+str(i) in chromTotRNAs: 
                                myChrom = myChrom+'.'+str(i)
                    myStart = int((splitLine[2]).split('=')[-1])
                    lineCounter = 0
                else:
                    myCoord = myChrom+'__'+str(myStart+lineCounter)
                    coordToPhyloP[myCoord] = float(splitLine[0])
                    lineCounter += 1

            tRNATotRNAPhyloPAvg = {}
            tRNATotRNAPhyloPAlign = {}
            tRNATo3PhyloPAvg = {}
            tRNATo3PhyloPAlign = {}
            tRNATo5PhyloPAvg = {}
            tRNATo5PhyloPAlign = {}

            for tRNA in tRNAToCoords:
                my5Val = 0.0
                my5Count = 0.0
                my3Val = 0.0
                my3Count = 0.0
                mytRNAVal = 0.0
                mytRNACount = 0.0
                for coord in tRNAToCoords[tRNA]:
                    mySprinzl = coordToSprinzl[coord]
                    if coord in coordToPhyloP:
                        if mySprinzl == 5:
                            my5Count += 1
                            my5Val += coordToPhyloP[coord]

                        elif mySprinzl == 3:
                            my3Count += 1
                            my3Val += coordToPhyloP[coord]

                        elif mySprinzl == 1:
                            mytRNACount += 1
                            mytRNAVal += coordToPhyloP[coord]

                if mytRNACount > 0:
                    tRNATotRNAPhyloPAvg[tRNA] = mytRNAVal/mytRNACount
                tRNATotRNAPhyloPAlign[tRNA] = mytRNACount
                if my5Count > 0:
                    tRNATo5PhyloPAvg[tRNA] = my5Val/my5Count
                tRNATo5PhyloPAlign[tRNA] = my5Count
                if my3Count > 0:
                    tRNATo3PhyloPAvg[tRNA] = my3Val/my3Count
                tRNATo3PhyloPAlign[tRNA] = my3Count

        ######################
        ######################
        ### GET BIT SCORES ###
        ######################
        ######################

        sys.stderr.write("Getting bit score features...\n")

        tRNAToGenBit = {}
        for line in open(self.inputOut):
            splitLine = (line.strip()).split()
            if not (splitLine[0]).strip() in ['Sequence', 'Name', '--------']:
                if splitLine[0]+'__'+str(min(int(splitLine[2]),int(splitLine[3]))+25) in coordTotRNA:
                    mytRNA = coordTotRNA[splitLine[0]+'__'+str(min(int(splitLine[2]),int(splitLine[3]))+25)]
                    if not mytRNA in myHaltRNAs:
                        mytRNA = myNameMap[mytRNA]
                    tRNAToGenBit[mytRNA] = float(splitLine[8])

        ###############
        ###############
        ### GET MFE ###
        ###############
        ###############

        tRNAToMFE = {}
        for line in open(self.inputMFE):
            splitLine = line.strip().split()
            if splitLine[0].startswith('>'):
                mytRNA = (splitLine[0])[1:]
            elif '(' in splitLine[0]:
                tRNAToMFE[mytRNA] = convertToFloat(splitLine[-1])

        ########################
        ########################
        ### GET CODON COUNTS ###
        ########################
        ########################

        tRNAToCodon = {}
        for tRNA in myHaltRNAs:
            myAC = tRNA[-3:]
            tRNAToCodon[tRNA] = ACtoCounts[myAC]

        ########################
        ########################
        ### GET CpG FEATURES ###
        ########################
        ########################

        sys.stderr.write("Getting sequence features...\n")

        tRNAToFasta = {}
        tRNAToCpGOvrPct = {}
        tRNAToObsExp = {}
        tRNAToObsExpUp = {}
        tRNAToObsExpDown = {}
        tRNAToTTTT = {}

        tRNAToUpstreamDist = {}
        tRNAToDownstreamDist = {}

        for line in open(self.inputBedExt):
            splitLine = (line.strip()).split('\t')
            myChrom = splitLine[0]
            mytRNA = splitLine[3]
            myStart = int(tRNAToStart[mytRNA].split('__')[-1])
            myEnd = int(tRNAToEnd[mytRNA].split('__')[-1])
            if splitLine[5] == '+':
                tRNAToUpstreamDist[mytRNA] = myStart-int(splitLine[1])
                tRNAToDownstreamDist[mytRNA] = int(splitLine[2])-myEnd
            else:
                tRNAToDownstreamDist[mytRNA] = myStart-int(splitLine[1])
                tRNAToUpstreamDist[mytRNA] = int(splitLine[2])-myEnd

        for line in open(self.inputFa):
            splitLine = (line.strip()).split('\t')
            mytRNA = splitLine[0]
            myFasta = splitLine[1]
            myUp = myFasta[:tRNAToUpstreamDist[mytRNA]]
            myDown = myFasta[-tRNAToDownstreamDist[mytRNA]:]
            tRNAToCpGOvrPct[mytRNA] = myFasta.count('CG')/float(len(myFasta))
            if (float(myFasta.count('C'))*float(myFasta.count('G'))) > 0:
                tRNAToObsExp[mytRNA] = (myFasta.count('CG')*float(len(myFasta))) / (float(myFasta.count('C'))*float(myFasta.count('G')))
            else:
                tRNAToObsExp[mytRNA] = 0.0

            if (float(myUp.count('C'))*float(myUp.count('G'))) > 0:
                tRNAToObsExpUp[mytRNA] = (myUp.count('CG')*float(len(myUp))) / (float(myUp.count('C'))*float(myUp.count('G')))
            else:
                tRNAToObsExpUp[mytRNA] = 0.0

            if (float(myDown.count('C'))*float(myDown.count('G'))) > 0:
                tRNAToObsExpDown[mytRNA] = (myDown.count('CG')*float(len(myDown))) / (float(myDown.count('C'))*float(myDown.count('G')))
            else:
                tRNAToObsExpDown[mytRNA] = 0.0
            
            if 'TTTT' in myDown:
                tRNAToTTTT[mytRNA] = myDown.index('TTTT')
            else:
                tRNAToTTTT[mytRNA] = len(myDown)+1


        ##################################
        ##################################
        ### GETTING PROXIMITY FEATURES ###
        ##################################
        ##################################

        sys.stderr.write("Getting proximity features...\n")

        chromToExonCoords = {}
        coordToExon = {}
        chromTotRNACoords = {}
        for chrom in chromTotRNAs:
            chromTotRNACoords[chrom] = []

        tRNATotRNA10kb = {}
        tRNAToProt75kb = {}

        if self.simplified == False:
            for line in open(self.inputGFF):
                splitLine = (line.strip()).split('\t')
                myChrom = splitLine[0]
                myExon = str(splitLine[3])
                if not myChrom in chromToExonCoords:
                    chromToExonCoords[myChrom] = [int(splitLine[1]), int(splitLine[2])]
                else:
                    (chromToExonCoords[myChrom]).append(int(splitLine[1]))
                    (chromToExonCoords[myChrom]).append(int(splitLine[2]))
                coordToExon[myChrom+'__'+str(splitLine[1])] = myExon
                coordToExon[myChrom+'__'+str(splitLine[2])] = myExon

            for myChrom in chromTotRNAs:
                for mytRNA in chromTotRNAs[myChrom]:
                    if myChrom in chromToExonCoords:
                        myExons75kb = set()
                        myStart = int(tRNAToStart[mytRNA].split('__')[1])
                        myEnd = int(tRNAToEnd[mytRNA].split('__')[1])
                        (chromTotRNACoords[myChrom]).append(myStart)
                        (chromTotRNACoords[myChrom]).append(myEnd)
                        for exonCoord in chromToExonCoords[myChrom]:
                            if min(abs(myStart-exonCoord), abs(myEnd-exonCoord)) > 10:
                                if min(abs(myStart-exonCoord), abs(myEnd-exonCoord)) <= 75000:
                                    myExons75kb.add(coordToExon[myChrom+'__'+str(exonCoord)])
                        tRNAToProt75kb[mytRNA] = len(myExons75kb)
                    else:
                        tRNAToProt75kb[mytRNA] = 0

        else:
            for myChrom in chromTotRNAs:
                for mytRNA in chromTotRNAs[myChrom]:
                    myStart = int(tRNAToStart[mytRNA].split('__')[1])
                    myEnd = int(tRNAToEnd[mytRNA].split('__')[1])
                    (chromTotRNACoords[myChrom]).append(myStart)
                    (chromTotRNACoords[myChrom]).append(myEnd)


        for myChrom in chromTotRNAs:
            for mytRNA in chromTotRNAs[myChrom]:
                myStart = int(tRNAToStart[mytRNA].split('__')[1])
                myEnd = int(tRNAToEnd[mytRNA].split('__')[1])
                mytRNAs10kb = set()
                for tRNACoord in chromTotRNACoords[myChrom]:
                    if not coordTotRNA[myChrom+'__'+str(tRNACoord)] == mytRNA:
                        addtRNA = coordTotRNA[myChrom+'__'+str(tRNACoord)]
                        if min(abs(myStart-tRNACoord), abs(myEnd-tRNACoord)) <= 10000:
                            mytRNAs10kb.add(addtRNA)
                tRNATotRNA10kb[mytRNA] = len(mytRNAs10kb)

        ############################
        ############################
        ### WRITE OUTPUT TO FILE ###
        ############################
        ############################

        if self.simplified == False:
            reducedSet = ['tRNAPhyloPAvg','5PhyloPAvg']
            reducedSet += ['CpGDensity','ObsExpUp','GenBit','tRNA10kb']
            reducedSet += ['Prot75kb','TTTT','Codon','MFE']
            reducedSetDicts = [tRNATotRNAPhyloPAvg,tRNATo5PhyloPAvg,tRNAToCpGOvrPct]
            reducedSetDicts += [tRNAToObsExpUp,tRNAToGenBit,tRNATotRNA10kb,tRNAToProt75kb,tRNAToTTTT,tRNAToCodon,tRNAToMFE]
        else:
            reducedSet = ['CpGDensity','ObsExpUp','GenBit','tRNA10kb','TTTT','Codon','MFE']
            reducedSetDicts = [tRNAToCpGOvrPct,tRNAToObsExpUp,tRNAToGenBit,tRNATotRNA10kb,tRNAToTTTT,tRNAToCodon,tRNAToMFE]

        tRNAToActivity = {}
        if len(self.activityFile) > 0:
            for line in open(self.activityFile):
                splitLine = (line.strip()).split()
                if str(splitLine[1]) in ['1', 'active']:
                    tRNAToActivity[splitLine[0]] = '1'
                elif str(splitLine[1]) in ['0', 'inactive']:
                    tRNAToActivity[splitLine[0]] = '0'
        else:
            for tRNA in myHaltRNAs:
                tRNAToActivity[tRNA] = '0'

        tRNAToSegDups = {}
        if len(self.segDups) > 0:
            for line in open(self.segDups):
                tRNAToSegDups[line.strip()] = True

        myOutString = 'tRNA\t'+joiner(reducedSet)+'\n'
        for tRNA in myHaltRNAs:
            if tRNA in tRNAToActivity:
                if not tRNA in tRNAToSegDups:
                    myOutString += tRNA
                    for diction in reducedSetDicts:
                        if not tRNA in diction:
                            diction[tRNA] = '?'
                        myOutString += '\t'+str(diction[tRNA])
                    myOutString += '\t'+tRNAToActivity[tRNA]+'\n'

        if self.simplified == False:
            myOutDataFile = 'alltRNAData.tsv'
        else:
            myOutDataFile = 'alltRNADataSimplified.tsv'
        open(myOutDataFile, 'w').write(myOutString)


        ###########################
        ###########################
        ### CLASSIFY TRNA GENES ###
        ###########################
        ###########################

        myHumanData = []
        myLabels = []
        myHumanNames = []
        imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
        for line in open(self.trainingFile):
            splitLine = (line.strip()).split('\t')
            if (splitLine[0]) != 'tRNA':
                myHumanData.append(makeFloat(splitLine[1:-1]))
                myHumanNames.append(splitLine[0])
                if str(splitLine[-1]) in ['active', '1']:
                    myLabels.append(1)
                elif str(splitLine[-1]) in ['inactive', '0']:
                    myLabels.append(0)

        imp_mean.fit_transform(myHumanData)
        myHumanDataReplaced = imp_mean.transform(myHumanData)

        clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
        clf.fit(myHumanDataReplaced, myLabels)

        myTestData = []
        myTestLabels = []
        myTestNames = []
        imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
        for line in open(myOutDataFile):
            splitLine = (line.strip()).split('\t')
            if (splitLine[0]) != 'tRNA':
                myTestData.append(makeFloat(splitLine[1:-1]))
                myTestNames.append(splitLine[0])

        imp_mean.fit_transform(myTestData)
        myTestDataReplaced = imp_mean.transform(myTestData)
        myScores = clf.predict_proba(myTestDataReplaced)
        myOutString = ''
        for i in range(0,len(myScores)):
            myOutString += myTestNames[i]+'\t'+getProb(myScores[i])+'\n'
        open(self.outputFile, 'w').write(myOutString)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ','.join(newList)

def convertToFloat(entry):
    myReturn = ''
    myNums = ['-','.','0','1','2','3','4','5','6','7','8','9']
    for i in range(0,len(entry)):
        if str(entry)[i] in myNums:
            myReturn += str(entry)[i]
    return(float(myReturn))

def makeFloat(myList):
    myReturn = []
    for k in myList:
        if not k == '?' and not 'tRNA' in k:
            myReturn.append(float(k))
        elif k == '?':
            myReturn.append(np.nan)
        else:
            myReturn.append(k)
    return(myReturn)

def getProb(myScores):
    if float(myScores[0]) > float(myScores[1]):
        return('-'+str(myScores[0])+'\t'+'inactive')
    else:
        return(str(myScores[1])+'\t'+'active')

def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    # Necessary files:
    if myCommandLine.args.inputBed:
        inputBed = myCommandLine.args.inputBed
    if myCommandLine.args.inputBedExt:
        inputBedExt = myCommandLine.args.inputBedExt
    if myCommandLine.args.inputOut:
        inputOut = myCommandLine.args.inputOut
    if myCommandLine.args.inputFa:
        inputFa = myCommandLine.args.inputFa
    if myCommandLine.args.inputMFE:
        inputMFE = myCommandLine.args.inputMFE
    if myCommandLine.args.trainingFile:
        trainingFile = myCommandLine.args.trainingFile

    # Files that we can deal with not having:
    if myCommandLine.args.phastConsElements:
        phastConsElements = myCommandLine.args.phastConsElements
    else:
        phastConsElements = ''

    if myCommandLine.args.inputWig:
        inputWig = myCommandLine.args.inputWig
    else:
        inputWig = ''

    if myCommandLine.args.inputGFF:
        inputGFF = myCommandLine.args.inputGFF
    else:
        inputGFF = ''

    if myCommandLine.args.chromLengths:
        chromLengths = myCommandLine.args.chromLengths
    else:
        chromLengths = ''

    if myCommandLine.args.activityFile:
        activityFile = myCommandLine.args.activityFile
    else:
        activityFile = ''

    if myCommandLine.args.segDups:
        segDups = myCommandLine.args.segDups
    else:
        segDups = ''

    # Flag for simplified version (no annotation/aligment) or not:
    if myCommandLine.args.simplified:
        simplified = True
    else:
        simplified = False

    # Specify output file
    # If nothing input by user, take .out tRNA file and add "tRNAPredictions" to that
    if myCommandLine.args.outputFile:
        outputFile = myCommandLine.args.outputFile
    if len(outputFile) == 0:
        outputFile = inputOut.split('.')[0]+'tRNAPredictions.out'

    myFileConverter = fileConverter(inputBed, inputBedExt, phastConsElements, inputWig, inputOut, inputMFE, inputFa, inputGFF, chromLengths, trainingFile, activityFile, segDups, outputFile, simplified)
    myFileConverter.getFeatures()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
