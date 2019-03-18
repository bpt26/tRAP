#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 6/6/2017
# gffToBed.py

import sys
import os
import time
import random
import numpy
import gzip
import math

"""
This program converts a .gff or .gff.gz file into a .bed file.
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
        self.parser.add_argument("-b", "--inputBed", help="Input .bed"+
            " file to be converted to a fasta.")
        self.parser.add_argument("-g", "--inputGenomeSeq", help="Input .fa containing "+
            "genome sequence to be parsed for sequences corresponding to .bed file.")
        self.parser.add_argument("-o", "--outputFile", help="The path to"+
            " and the filename of the .fa file you are creating", default='')
        self.args = self.parser.parse_args()


class fileConverter(object):
    """
    Primary class where filetype is converted.
    """

    def __init__(self, inputBed, inputGenomeSeq, outputFile):
        self.inputBed = inputBed
        self.inputGenomeSeq = inputGenomeSeq
        self.outputFile = outputFile


    def convertFile(self):
        chromTotRNACoords = {}
        for line in open(self.inputBed):
            splitLine = (line.strip()).split('\t')
            myChrom = splitLine[0]
            myStart = int(splitLine[1])
            myEnd = int(splitLine[2])
            myName = str(splitLine[3])
            myStrand = str(splitLine[5])
            if not myChrom in chromTotRNACoords:
                chromTotRNACoords[myChrom] = []
            (chromTotRNACoords[myChrom]).append([myStart, myEnd, myName, myStrand])

        currentChrom = ''
        myCurrentIndex = 0
        tRNAToSeq = {}
        for line in open(self.inputGenomeSeq):
            stripLine = line.strip()
            if stripLine.startswith('>'):
                if len(currentChrom) > 0:
                    for tRNA in mytRNAs:
                        if tRNA[3] == '+':
                            tRNAToSeq[str(tRNA[2])] = currentChrom[(tRNA[0]):(tRNA[1])]
                        else:
                            tRNAToSeq[str(tRNA[2])] = revComp(currentChrom[(tRNA[0]):(tRNA[1])])

                myChrom = stripLine[1:]
                currentChrom = ''
                mytRNAs = []
                if myChrom in chromTotRNACoords:
                    mytRNAs = sorted(chromTotRNACoords[myChrom], key=lambda x: x[1])
            elif len(mytRNAs) > 0:
                currentChrom += (stripLine.upper())

        if len(currentChrom) > 0:
            for tRNA in mytRNAs:
                if tRNA[3] == '+':
                    tRNAToSeq[str(tRNA[2])] = currentChrom[(tRNA[0]):(tRNA[1])]
                else:
                    tRNAToSeq[str(tRNA[2])] = revComp(currentChrom[(tRNA[0]):(tRNA[1])])

        seqTotRNA = {}
        myOutString = ''
        for tRNA in sorted(tRNAToSeq.keys()):
            mySeq = tRNAToSeq[tRNA]
            myOutString += tRNA+'\t'+mySeq+'\n'
            if not mySeq in seqTotRNA:
                seqTotRNA[mySeq] = tRNA
        open(self.outputFile, 'w').write(myOutString)

def revComp(seq):
    diction = {}
    diction['A'] = 'T'
    diction['C'] = 'G'
    diction['G'] = 'C'
    diction['T'] = 'A'
    diction['N'] = 'N'
    mySeq = seq[::-1]
    myReturn = ''
    i = 0
    while i < len(mySeq):
        myReturn += diction[mySeq[i]]
        i += 1
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)


def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    if myCommandLine.args.inputFile:
        inputFile = myCommandLine.args.inputFile

    if myCommandLine.args.outputFile:
        outputFile = myCommandLine.args.outputFile

    if len(myCommandLine.args.outputFile) == 0:
        outputFile = inputFile.split('.')[0]+'.bed'

    myFileConverter = fileConverter(inputFile, outputFile)
    myFileConverter.convertFile()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit