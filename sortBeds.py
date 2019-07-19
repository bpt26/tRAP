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

def sortBeds():
	myOutString = ''
	startToLine = {}
	for line in open(str(sys.argv[1])):
		splitLine = (line.strip()).split('\t')
		startToLine[int(splitLine[1])] = splitLine
	for k in sorted(startToLine.keys()):
		myOutString += joiner(startToLine[k])+'\n'
	open(str(sys.argv[1]),'w').write(myOutString)

def joiner(entry):
    """
    Helper function to print lists in 
    tab separated format.
    """
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def main():
    sortBeds()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

