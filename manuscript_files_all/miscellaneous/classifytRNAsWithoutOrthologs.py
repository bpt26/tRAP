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


def classifyNoOrthologs():
    mouseNonHuman = {}
    counter = 0
    for line in open('FINAL_ORTHOLOG_SET.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[13] == 'N/A' and splitLine[18] != 'N/A':
            mouseNonHuman[(splitLine[18]).split('Mmusc10-')[-1]] = True
        elif splitLine[13] != 'N/A' and splitLine[18] != 'N/A':
            counter += 1
    print(len(mouseNonHuman))
    print(counter)


    myHumanData = []
    myLabels = []
    myHumanNames = []
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    for line in open('humanCpGTrainingSetFixed.tsv'):
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
    for line in open('mousetRNADataCpGFixed.tsv'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] in mouseNonHuman:
            myTestData.append(makeFloat(splitLine[1:-1]))
            myTestNames.append(splitLine[0])
            myTestLabels.append(str(splitLine[-1]))

    imp_mean.fit_transform(myTestData)
    myTestDataReplaced = imp_mean.transform(myTestData)
    myScores = clf.predict_proba(myTestDataReplaced)

    myPredictions = clf.predict(myTestDataReplaced)
    for i in range(0,len(myPredictions)):
        if not str(myPredictions[i]) == str(myTestLabels[i]):
            print(myTestNames[i], str(myPredictions[i]), str(myTestLabels[i]))
    cM = confusionMatrix(clf.predict(myTestDataReplaced),myTestLabels)
    print(cM)
    print(getAccuracy(cM))
    print(getScore(clf.predict_proba(myTestDataReplaced),myTestLabels))



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

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\n'.join(newList)

def getScore(list1,list2):
    myTotal = 0.0
    for i in range(0,len(list1)):
        if ((list1[i])[0] > (list1[i])[1] and int(list2[i]) == 0) or ((list1[i])[1] > (list1[i])[0] and int(list2[i]) == 1):
            myTotal += 1.0
    return(myTotal/float(len(list2)))

def confusionMatrix(pred,real):
    myReturn = [[0.0, 0.0],[0.0, 0.0]]
    for i in range(0,len(pred)):
        if int(pred[i]) == 0 and int(real[i]) == 0:
            (myReturn[0])[0] += 1
        elif int(pred[i]) == 0 and int(real[i]) == 1:
            (myReturn[0])[1] += 1
        elif int(pred[i]) == 1 and int(real[i]) == 0:
            (myReturn[1])[0] += 1
        elif int(pred[i]) == 1 and int(real[i]) == 1:
            (myReturn[1])[1] += 1
    return(myReturn)

def getAccuracy(cM):
    return( (((cM[0])[0]+(cM[1])[1]) / ((cM[0])[0]+(cM[0])[1]+(cM[1])[0]+(cM[1])[1])) * 100.0 )

def getFirst(myList):
    myReturn = []
    for k in myList:
        myReturn.append(k[0])
    return(myReturn)



def main():
    classifyNoOrthologs()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit




