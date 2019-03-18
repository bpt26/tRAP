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
        self.parser.add_argument("-i", "--inputFile", help="Input .gff"+
            " or .gff.gz file.")
        self.parser.add_argument("-o", "--outputFile", help="The path to"+
            " and the filename of the bed file you are creating", default='')
        self.args = self.parser.parse_args()


class fileConverter(object):
    """
    Primary class where filetype is converted.
    """

    def __init__(self, inputFile, outputFile):
        self.inputFile = inputFile
        self.outputFile = outputFile

    def convertFile(self):
        myOutString = ''
        open(self.outputFile, 'w').write('')
        alreadyPrinted = {}
        if self.inputFile.endswith('gz'):
            for line in gzip.open(self.inputFile):
                newLine = self.convertLine((line.strip()).split('\t'))
                if len(newLine) > 2 and not joiner(newLine[:2]) in alreadyPrinted and not joiner([newLine[0],newLine[2]]) in alreadyPrinted:
                    myOutString += joiner(newLine)+'\n'
                    alreadyPrinted[joiner(newLine[:2])] = True
                    alreadyPrinted[joiner([newLine[0],newLine[2]])] = True
                    if len(alreadyPrinted.keys()) % 10000 == 0:
                        sys.stderr.write(str(len(alreadyPrinted.keys()))+' lines processed!\n')
                        open(self.outputFile, 'a').write(myOutString)
                        myOutString = ''
        else:
            for line in open(self.inputFile):
                newLine = self.convertLine((line.strip()).split('\t'))
                if len(newLine) > 2 and not joiner(newLine[:2]) in alreadyPrinted and not joiner([newLine[0],newLine[2]]) in alreadyPrinted:
                    myOutString += joiner(newLine)+'\n'
                    alreadyPrinted[joiner(newLine[:2])] = True
                    alreadyPrinted[joiner([newLine[0],newLine[2]])] = True
                    if len(alreadyPrinted.keys()) % 10000 == 0:
                        sys.stderr.write(str(len(alreadyPrinted.keys()))+' lines processed!\n')
                        open(self.outputFile, 'a').write(myOutString)
                        myOutString = ''
        open(self.outputFile, 'a').write(myOutString)

        

    def convertLine(self, splitLine):
        # We want *only* exons, and tRNAs are sometimes labeled as "exons annotated by tRNAscan-SE"
        if (len(splitLine) > 2) and (str(splitLine[2]) == 'exon' and not 'trna' in str(splitLine[1]).lower()):
            # Adjust because GFFs are 1-based and BEDs are 0-based:
            myStart = int(splitLine[3])-1
            myEnd = int(splitLine[4])-1
            myInfo = splitLine[8].split(';')

            # Cover differences in format:
            if ':' in str(myInfo[2]):
                myGene = str((myInfo[2]).split(':')[-1])
            else:
                myGene = str((myInfo[2]).split('=')[-1])
            myGene = ''.join(myGene.split(' '))

            if ':' in str(myInfo[0]):
                myID = str((myInfo[0]).split(':')[-1])
            else:
                myID = str((myInfo[0]).split('=')[-1])
            myID = ''.join(myID.split(' '))

            #print([splitLine[0],myStart,myEnd,myGene+'.'+myID,0,splitLine[6],myStart,myEnd,'255,0,0\t1',str(myEnd-myStart),0])
            return([splitLine[0],myStart,myEnd,myGene+'.'+myID,0,splitLine[6],myStart,myEnd,'255,0,0\t1',str(myEnd-myStart),0])
        else:
            return([])

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

