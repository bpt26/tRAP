#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 6/6/2017
# makeHiConfBeds.py

import sys
import os
import time
import random
import numpy
import gzip
import math

"""
This file outputs three main things:
    1) the bed file corresponding to all high-confidence tRNA genes
    2) a folder of bed files containing the high-confidence tRNA genes separated by chromosome
    3) a .sh script that will make .mafs out of all files in 2) and all of the .wig and .elements files
"""

import sys, argparse, random, math

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
        self.parser.add_argument("-i", "--hiConfOut", help="Input .out"+
            " of high-confidence tRNAs from tRNAscan-SE 2.0.")
        self.parser.add_argument("-b", "--bedFile", help="Input .bed"+
            " of all tRNAs from tRNAscan-SE 2.0.")
        self.parser.add_argument("-s", "--speciesName", help="Name of species"+
            " of interest as it appears in the Cactus graph.", default='')
        self.parser.add_argument("-c", "--cactusPath", help="Path to Cactus "+
            " graph or .hal file.", default='')
        self.parser.add_argument("-m", "--modFile", help="Path to .mod "+
            " file that will be used to train the PhyloP and PhastCons models.", default='')
        self.parser.add_argument("-l", "--chromLengths", help="File containing"+
            " length of each chromosome corresponding to this genome.", default='')
        self.parser.add_argument("-o", "--outPath", help="Path to your desired"+
            " output including prefix for file names.", default='tRNAHiConf')
        self.args = self.parser.parse_args()


class fileConverter(object):
    """
    Primary class where filetype is converted.
    """
    def __init__(self, hiConfOut, bedFile, speciesName, cactusPath, modFile, chromLengths, outPath):
        self.hiConfOut = hiConfOut
        self.bedFile = bedFile
        self.speciesName = speciesName
        self.cactusPath = cactusPath
        self.modFile = modFile
        self.chromLengths = chromLengths
        self.outPath = outPath

    def makeHiConfBeds(self):
        coordTotRNA = {}
        for line in open(self.hiConfOut):
            splitLine = (line.strip()).split()
            if not (splitLine[0].strip()) in ['Sequence','Name','--------']:
                for k in range(min(int(splitLine[2]),int(splitLine[3])),max(int(splitLine[2]),int(splitLine[3]))+1):
                    coordTotRNA[str(splitLine[0])+'-'+str(k)] = True

        chromToLengths = {}
        if len(self.chromLengths) > 0:
            for line in open(self.chromLengths):
                splitLine = (line.strip()).split('\t')
                if len(splitLine) > 2:
                    chromToLengths[str(splitLine[0])] = int(splitLine[2])

        myOutString = ''
        myOutString50 = ''
        myOutString250 = ''
        chromToOut = {}

        for line in open(self.bedFile):
            splitLine = (line.strip()).split('\t')
            if splitLine[0]+'-'+str(int(splitLine[1])+15) in coordTotRNA:

                myOutString += joiner(splitLine)+'\n'

                if not splitLine[0] in chromToOut:
                    chromToOut[splitLine[0]] = []

                if splitLine[0] in chromToLengths:
                    (chromToOut[splitLine[0]]).append(makeNewLineAdd(splitLine,chromToLengths[splitLine[0]],50)+'\n')
                    myOutString50 += makeNewLineAdd(splitLine,chromToLengths[splitLine[0]],50)+'\n'
                    myOutString250 += makeNewLineAdd(splitLine,chromToLengths[splitLine[0]],250)+'\n'
                else:
                    (chromToOut[splitLine[0]]).append(makeNewLineAdd(splitLine,99999999999999,50)+'\n')
                    myOutString50 += makeNewLineAdd(splitLine,99999999999999,50)+'\n'
                    myOutString250 += makeNewLineAdd(splitLine,99999999999999,250)+'\n'

        open(self.outPath+'.bed', 'w').write(myOutString)
        open(self.outPath+'50.bed', 'w').write(myOutString50)
        open(self.outPath+'250.bed', 'w').write(myOutString250)

        myCommands = '#!/bin/bash\n#$ -cwd\n#$ -j y\n#$ -S /bin/bash\n\nmkdir chromWigs\nmkdir chromElements\n'
        for chrom in sorted(chromToOut.keys()):
            # unfortunately haven't found a way around just running as 
            myCommands += 'rm *temp*\nrm *Temp*\nhal2maf --refGenome '+self.speciesName
            myCommands += ' --refTargets '+chrom+'.bed '+self.cactusPath+' '+chrom+'.maf\n'
            myCommands += 'phyloP --msa-format MAF --method LRT --wig-scores --mode CONACC --no-prune '+self.modFile
            myCommands += ' '+chrom+'.maf > chromWigs/'+chrom+'.wig\nphastCons '+chrom+'.maf '+self.modFile+' --msa-format MAF --viterbi '
            myCommands += 'chromElements/'+chrom+'_phastCons.bed --score\n\n'
            open(chrom+'.bed', 'w').write(''.join(chromToOut[chrom]))
        open('getAlign.sh', 'w').write(myCommands+'cat chromWigs/*.wig > '+self.speciesName+'.wig\ncat chromElements/*phastCons.bed > ')
        open('getAlign.sh', 'a').write(self.speciesName+'_elements.bed\n')


def joiner(entry):
    """
    Helper function to print lists in 
    tab separated format.
    """
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)


def makeNewLineAdd(oldLine, myLen, distance):
    """
    Adds distance on either side of line of bed, making sure
    not to go over length of chromosome.
    """

    oldLine[1] = int(oldLine[1])
    oldLine[2] = int(oldLine[2])
    oldLine[6] = int(oldLine[6])
    oldLine[7] = int(oldLine[7])

    if oldLine[1] <= int(distance):
        oldLine[1] = 0
        oldLine[6] = 0
    else:
        oldLine[1] -= distance
        oldLine[6] -= distance

    if oldLine[2]+distance >= myLen:
        oldLine[2] = myLen-1
        oldLine[7] = myLen-1
    else:
        oldLine[2] += distance
        oldLine[7] += distance

    oldLine[9] = '1'
    oldLine[10] = str(oldLine[2]-oldLine[1])+','
    oldLine[11] = '0,'
    return(joiner(oldLine))


def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    if myCommandLine.args.hiConfOut:
        hiConfOut = myCommandLine.args.hiConfOut

    if myCommandLine.args.bedFile:
        bedFile = myCommandLine.args.bedFile

    if myCommandLine.args.speciesName:
        speciesName = myCommandLine.args.speciesName
    else:
        speciesName = ''

    if myCommandLine.args.cactusPath:
        cactusPath = myCommandLine.args.cactusPath
    else:
        cactusPath = ''

    if myCommandLine.args.modFile:
        modFile = myCommandLine.args.modFile
    else:
        modFile = ''

    if myCommandLine.args.chromLengths:
        chromLengths = myCommandLine.args.chromLengths
    else:
        chromLengths = ''

    if myCommandLine.args.outPath:
        outPath = myCommandLine.args.outPath

    myFileConverter = fileConverter(hiConfOut, bedFile, speciesName, cactusPath, modFile, chromLengths, outPath)
    myFileConverter.makeHiConfBeds()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
