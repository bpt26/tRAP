#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 10/16/2018
# newerParseBlast.py

import sys
import os
import time
import random
import numpy
import gzip
import math

##########################
##### MAIN FUNCTIONS #####
##########################

def get95CI():
    catToActiveVec = {}
    catToInactiveVec = {}
    for line in open('humanCpGTrainingSet.tsv'):
        splitLine = (line.strip()).split('\t')
        if str(splitLine[-1]) == '1':
            for i in range(1,len(splitLine)-1):
                if not splitLine[i] == '?':
                    catToActiveVec[myHeader[i]].append(float(splitLine[i]))
        elif str(splitLine[-1]) == '0':
            for i in range(1,len(splitLine)-1):
                if not splitLine[i] == '?':
                    catToInactiveVec[myHeader[i]].append(float(splitLine[i]))
        elif splitLine[-1] == 'activity':
            myHeader = splitLine
            for i in range(1,len(splitLine)-1):
                catToActiveVec[splitLine[i]] = []
                catToInactiveVec[splitLine[i]] = []

    # print(catToActiveVec)
    # print(catToInactiveVec)

    for catNum in range(1,len(myHeader)-1):
        # ACTIVE:
        myMeans = []
        for iterNum in range(0,1000):
            myTempList = []
            for i in range(0,len(catToActiveVec[myHeader[catNum]])):
                myTempList.append(random.sample(catToActiveVec[myHeader[catNum]], 1)[0])
            myMeans.append(getMean(myTempList))
        myMeans.sort()
        print(myHeader[catNum], 'ACTIVE', getMean(catToActiveVec[myHeader[catNum]]), myMeans[24], myMeans[974])
        # INACTIVE:
        myMeans = []
        for iterNum in range(0,1000):
            myTempList = []
            for i in range(0,len(catToInactiveVec[myHeader[catNum]])):
                myTempList.append(random.sample(catToInactiveVec[myHeader[catNum]], 1)[0])
            myMeans.append(getMean(myTempList))
        myMeans.sort()
        print(myHeader[catNum], 'INACTIVE', getMean(catToInactiveVec[myHeader[catNum]]), myMeans[24], myMeans[974])

def getMean(myList):
    myReturn = 0.0
    for k in myList:
        myReturn += float(k)
    return(myReturn/float(len(myList)))





def main():
    get95CI()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit


