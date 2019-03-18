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
from random import shuffle
import operator

"""
This program grabs a set number of exons from a bed file of 4d sites
to make a smaller file more easily used to train a PhyloP model.
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
        self.parser.add_argument("-i", "--inputFile", help="Input .bed"+
            " of 4d sites to reduce.")
        self.parser.add_argument("-o", "--outputFile", help="The path to"+
            " your desired output file.", default='')
        self.parser.add_argument("-n", "--numExons", help="Number of"+
            " exons to include in the file.", default=100000)
        self.args = self.parser.parse_args()

class fileConverter(object):
    """
    Primary class where filetype is converted.
    """

    def __init__(self, inputFile, outputFile, numExons):
        self.inputFile = inputFile
        self.outputFile = outputFile
        self.numExons = numExons

    def reduce4dBeds(self):
        allLines = []
        keepLines = []
        myOutString = ''
        for line in open(self.inputFile):
            allLines.append((line.strip()).split('\t'))
        shuffle(allLines)
        for i in range(0,self.numExons):
            keepLines.append(allLines[i])
        for k in sorted(keepLines, key=operator.itemgetter(0,1)):
            myOutString += joiner(k)+'\n'
        open(self.outputFile, 'w').write(myOutString)

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
        outputFile = inputFile.split('.')[0]+'reduced.bed'

    if myCommandLine.args.numExons:
        numExons = int(myCommandLine.args.numExons)

    myFileConverter = fileConverter(inputFile, outputFile, numExons)
    myFileConverter.reduce4dBeds()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit







