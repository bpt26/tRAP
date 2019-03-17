#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 6/6/2017
# getMFE.py

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
        self.parser.add_argument("-s", "--ssFile", help="Input .ss"+
            " of high-confidence tRNAs from tRNAscan-SE 2.0.")
        self.parser.add_argument("-b", "--bedFile", help="Input .bed"+
            " of high-confidence tRNAs from tRNAscan-SE 2.0.")
        self.parser.add_argument("-o", "--outFile", help="File to be analyzed"+
            " by RNAfold to get MFE data for tRNAs.", default='tRNAFoldIn.txt')
        self.args = self.parser.parse_args()


class fileConverter(object):
    """
    Primary class where filetype is converted.
    """

    def __init__(self, ssFile, bedFile, outFile):
        self.ssFile = ssFile
        self.bedFile = bedFile
        self.outFile = outFile

    def processStr(self, myString):
        toReturn = ''
        for i in range(0,len(myString)):
            if myString[i] == '>':
                toReturn += '('
            elif myString[i] == '.':
                toReturn += '.'
            elif myString[i] == '<':
                toReturn += ')'
        return(toReturn)

    def convertFile(self):
        myOutString = ''
        coordTotRNA = {}
        for line in open(self.bedFile):
            splitLine = (line.strip()).split('\t')
            for k in range(int(splitLine[1]),int(splitLine[2])+1):
                coordTotRNA[splitLine[0]+'__'+str(k)] = str(splitLine[3])
        mytRNA = ''
        for line in open(self.ssFile):
            splitLine = (line.strip()).split()
            # first line
            if len(splitLine) > 2 and 'Length' in splitLine[2]:
                myChrom = splitLine[0].split('.trna')[0]
                myCoords = ((splitLine[1])[1:-1]).split('-')
                myStart = min(int(myCoords[0]), int(myCoords[1]))
                myEnd = min(int(myCoords[0]), int(myCoords[1]))
                #print(myStart, myEnd)
                if myChrom+'__'+str(myEnd+25) in coordTotRNA:
                    mytRNA = coordTotRNA[myChrom+'__'+str(myEnd+25)]
            # seq line
            if len(mytRNA) > 1 and (splitLine[0]).startswith('Seq:'):
                mySeq = (splitLine[1]).upper()
            if len(mytRNA) > 1 and (splitLine[0]).startswith('Str:'):
                myOutString += '>'+mytRNA+'\n'+mySeq+'\n'+self.processStr(splitLine[1])+'\n'
                mytRNA = ''
        open(self.outFile, 'w').write(myOutString)


def main(myCommandLine=None):
    """
    Initializes a CommandLine object and passes the provided 
    arguments into a new fileConverter object and calls main method.
    """
    myCommandLine = CommandLine()

    if myCommandLine.args.ssFile:
        ssFile = myCommandLine.args.ssFile

    if myCommandLine.args.bedFile:
        bedFile = myCommandLine.args.bedFile

    if myCommandLine.args.outFile:
        outFile = myCommandLine.args.outFile

    myFileConverter = fileConverter(ssFile, bedFile, outFile)
    myFileConverter.convertFile()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit