#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 10/16/2018
# oneToOneMap.py

import sys
import os
import time
import random
import numpy
import gzip
import math
import scipy
from scipy import stats

##########################
##### MAIN FUNCTIONS #####
##########################

def correlateFeatures():
    featureToList = {}
    for line in open('humanCpGTrainingSet.tsv'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'tRNA':
            myHeader = splitLine
            for i in range(1,len(splitLine)):
                featureToList[i] = []
        elif not '?' in splitLine:
            for i in range(1,len(splitLine)):
                featureToList[i].append(float(splitLine[i]))

    myOutString = ''
    for i in range(1,12):
        for m in range(1,12):
            if i != m:
                myOutString += myHeader[i]+'\t'+myHeader[m]+'\t'+str(scipy.stats.spearmanr(featureToList[i],featureToList[m])[0])+'\t'+str(scipy.stats.spearmanr(featureToList[i],featureToList[m])[1])+'\n'
    open('correlatedFeatures.txt', 'w').write(myOutString)

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
    correlateFeatures()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit