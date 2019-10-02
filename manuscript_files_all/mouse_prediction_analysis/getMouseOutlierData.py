#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 10/16/2018
# oneToOneMap.py

import sys
import os
import time
import random
import numpy as np
import gzip
import math
import scipy
from scipy import stats

##########################
##### MAIN FUNCTIONS #####
##########################

### NOTE: RUN PLOTMOUSEBYGROUP.PY BEFORE YOU DO THIS

def getMouseOutlierData():
    tRNAToConf = {}
    tRNAToGroup = {}
    tRNAToCount = {}
    tRNAToActivity = {}
    colors = ['red','darkorange','gold','green','blue']
    groups = ['A','B','C','D','E']

    for line in open('mouse_liver_1.bed'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            tRNAToConf[splitLine[0]] = float(splitLine[1])
            tRNAToGroup[splitLine[0]] = groups[colors.index(splitLine[2])]
            tRNAToCount[splitLine[0]] = float(splitLine[3])
            if tRNAToGroup[splitLine[0]] in ['A','B','C']:
                tRNAToActivity[splitLine[0]] = 'active'
            else:
                tRNAToActivity[splitLine[0]] = 'inactive'

    featureToActiveData = {}
    featureToInactiveData = {}

    for line in open('humanCpGTrainingSet.tsv'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'tRNA':
            myHeader = splitLine
        else:
            if int(splitLine[-1]) == 0:
                for i in range(1,len(splitLine)-1):
                    if not str(splitLine[i]) == '?':
                        if not myHeader[i] in featureToInactiveData:
                            featureToInactiveData[myHeader[i]] = []
                        featureToInactiveData[myHeader[i]].append(float(splitLine[i]))
            elif int(splitLine[-1]) == 1:
                for i in range(1,len(splitLine)-1):
                    if not str(splitLine[i]) == '?':
                        if not myHeader[i] in featureToActiveData:
                            featureToActiveData[myHeader[i]] = []
                        featureToActiveData[myHeader[i]].append(float(splitLine[i]))

    featureToActiveMean = {}
    featureToActiveSD = {}
    featureToInactiveMean = {}
    featureToInactiveSD = {}
    for feature in myHeader[1:-1]:
        featureToActiveMean[feature] = np.mean(featureToActiveData[feature])
        featureToActiveSD[feature] = np.std(featureToActiveData[feature])
        featureToInactiveMean[feature] = np.mean(featureToInactiveData[feature])
        featureToInactiveSD[feature] = np.std(featureToInactiveData[feature])
        print(feature, featureToActiveSD[feature], featureToInactiveSD[feature])
        print()

    # myOutliers = ['NC_000077.6.tRNA94-AspGTC','NC_000077.6.tRNA1235-CysGCA','NC_000079.6.tRNA90-AlaTGC']
    # myOutliers += ['NC_000079.6.tRNA1532-ThrAGT','NC_000079.6.tRNA1531-GlnTTG','NC_000079.6.tRNA1525-LeuTAA']
    # myOutliers += ['NC_000079.6.tRNA110-GlnCTG','NC_000079.6.tRNA111-SerAGA','NC_000079.6.tRNA1517-iMetCAT']
    # myOutliers += ['NC_000079.6.tRNA1515-AspGTC','NC_000079.6.tRNA1505-ThrAGT','NC_000079.6.tRNA122-MetCAT']
    # myOutliers += ['NC_000079.6.tRNA125-TyrGTA','NC_000079.6.tRNA129-GluTTC','NC_000079.6.tRNA132-iMetCAT']
    # myOutliers += ['NC_000069.6.tRNA463-GluTTC','NC_000069.6.tRNA1203-GlyCCC','NC_000069.6.tRNA464-GlyCCC']
    # myOutliers += ['NC_000069.6.tRNA1202-GluTTC','NC_000069.6.tRNA1201-GlnCTG','NC_000069.6.tRNA470-LysCTT']
    # myOutliers += ['NC_000069.6.tRNA472-AsnGTT','NC_000069.6.tRNA473-HisGTG','NC_000069.6.tRNA1199-HisGTG']
    # myOutliers += ['NC_000069.6.tRNA1198-AsnGTT','NC_000070.6.tRNA979-AlaAGC','NC_000086.7.tRNA1530-GlnTTG']
    # myOutliers += ['NC_000086.7.tRNA1043-LysCTT','NC_000086.7.tRNA604-AlaTGC','NC_000086.7.tRNA606-AlaTGC','NC_000086.7.tRNA613-AlaTGC']

    myOutliers = ['NC_000077.6.tRNA94-AspGTC','NC_000077.6.tRNA1235-CysGCA','NC_000079.6.tRNA90-AlaTGC']
    myOutliers += ['NC_000079.6.tRNA1532-ThrAGT','NC_000079.6.tRNA1531-GlnTTG','NC_000079.6.tRNA1525-LeuTAA']
    myOutliers += ['NC_000079.6.tRNA110-GlnCTG','NC_000079.6.tRNA111-SerAGA','NC_000079.6.tRNA1517-iMetCAT']
    myOutliers += ['NC_000079.6.tRNA1515-AspGTC','NC_000079.6.tRNA1505-ThrAGT','NC_000079.6.tRNA122-MetCAT']
    myOutliers += ['NC_000079.6.tRNA125-TyrGTA','NC_000079.6.tRNA129-GluTTC','NC_000079.6.tRNA132-iMetCAT']
    myOutliers += ['NC_000069.6.tRNA463-GluTTC','NC_000069.6.tRNA1203-GlyCCC','NC_000069.6.tRNA464-GlyCCC']
    myOutliers += ['NC_000069.6.tRNA1202-GluTTC','NC_000069.6.tRNA1201-GlnCTG','NC_000069.6.tRNA470-LysCTT']
    myOutliers += ['NC_000069.6.tRNA472-AsnGTT','NC_000069.6.tRNA473-HisGTG','NC_000069.6.tRNA1199-HisGTG']
    myOutliers += ['NC_000069.6.tRNA1198-AsnGTT','NC_000070.6.tRNA979-AlaAGC','NC_000086.7.tRNA1530-GlnTTG']
    myOutliers += ['NC_000086.7.tRNA1043-LysCTT','NC_000086.7.tRNA604-AlaTGC','NC_000086.7.tRNA606-AlaTGC','NC_000086.7.tRNA613-AlaTGC']

    outlierDict = {'NC_000077.6.tRNA94-AspGTC': 'tRNA-Asp-GTC-1-9',
        'NC_000077.6.tRNA1235-CysGCA': 'tRNA-Cys-GCA-3-2',
        'NC_000079.6.tRNA90-AlaTGC': 'tRNA-Ala-TGC-1-1',
        'NC_000079.6.tRNA1532-ThrAGT': 'tRNA-Thr-AGT-6-1',
        'NC_000079.6.tRNA1531-GlnTTG': 'tRNA-Gln-TTG-2-1',
        'NC_000079.6.tRNA1525-LeuTAA': 'tRNA-Leu-TAA-3-1',
        'NC_000079.6.tRNA110-GlnCTG': 'tRNA-Gln-CTG-1-1',
        'NC_000079.6.tRNA111-SerAGA': 'tRNA-Ser-AGA-2-4',
        'NC_000079.6.tRNA1517-iMetCAT': 'tRNA-iMet-CAT-1-3',
        'NC_000079.6.tRNA1515-AspGTC': 'tRNA-Asp-GTC-1-12',
        'NC_000079.6.tRNA1505-ThrAGT': 'tRNA-Thr-AGT-5-1',
        'NC_000079.6.tRNA122-MetCAT': 'tRNA-Met-CAT-4-1',
        'NC_000079.6.tRNA125-TyrGTA': 'tRNA-Tyr-GTA-1-2',
        'NC_000079.6.tRNA129-GluTTC': 'tRNA-Glu-TTC-1-3',
        'NC_000079.6.tRNA132-iMetCAT': 'tRNA-iMet-CAT-1-5',
        'NC_000069.6.tRNA463-GluTTC': 'tRNA-Glu-TTC-3-1',
        'NC_000069.6.tRNA1203-GlyCCC': 'tRNA-Gly-CCC-2-1',
        'NC_000069.6.tRNA464-GlyCCC': 'tRNA-Gly-CCC-2-2',
        'NC_000069.6.tRNA1202-GluTTC': 'tRNA-Glu-TTC-3-2',
        'NC_000069.6.tRNA1201-GlnCTG': 'tRNA-Gln-CTG-3-2',
        'NC_000069.6.tRNA470-LysCTT': 'tRNA-Lys-CTT-3-1',
        'NC_000069.6.tRNA472-AsnGTT': 'tRNA-Asn-GTT-3-4',
        'NC_000069.6.tRNA473-HisGTG': 'tRNA-His-GTG-2-5',
        'NC_000069.6.tRNA1199-HisGTG': 'tRNA-His-GTG-2-6',
        'NC_000069.6.tRNA1198-AsnGTT': 'tRNA-Asn-GTT-3-5',
        'NC_000070.6.tRNA979-AlaAGC': 'tRNA-Ala-AGC-8-1',
        'NC_000086.7.tRNA1530-GlnTTG': 'tRNA-Gln-TTG-4-1',
        'NC_000086.7.tRNA1043-LysCTT': 'tRNA-Lys-CTT-1-1',
        'NC_000086.7.tRNA604-AlaTGC': 'tRNA-Ala-TGC-5-1',
        'NC_000086.7.tRNA606-AlaTGC': 'tRNA-Ala-TGC-5-2',
        'NC_000086.7.tRNA613-AlaTGC': 'tRNA-Ala-TGC-6-1'}

    for line in open('mousetRNADataCpG.tsv'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'tRNA':
            myOutString = joiner(splitLine)+'\tProbability Score\tExperimental Data\tClassifier Prediction\n'
            activeMeans = []
            activeSDs = []
            inactiveMeans = []
            inactiveSDs = []
            myHeader = splitLine
            for feature in myHeader[1:]:
                activeMeans.append(featureToActiveMean[feature])
                activeSDs.append(featureToActiveSD[feature])
                inactiveMeans.append(featureToInactiveMean[feature])
                inactiveSDs.append(featureToInactiveSD[feature])
            myOutString += 'activeMeans\t'+joiner(activeMeans)+'\nactiveSDs\t'+joiner(activeSDs)+'\n'
            myOutString += 'inactiveMeans\t'+joiner(inactiveMeans)+'\ninactiveSDs\t'+joiner(inactiveSDs)+'\n\n'
        elif splitLine[0] in myOutliers:
            myOutString += outlierDict[splitLine[0]]+'\t'+joiner(splitLine[1:-1])+'\t'+str(tRNAToConf[splitLine[0]])+'\t'+tRNAToActivity[splitLine[0]]
            if tRNAToActivity[splitLine[0]] == 'active':
                myOutString += '\tinactive\n'
            else:
                myOutString += '\tactive\n'
            myActiveZScores = []
            for i in range(1,len(myHeader)):
                myActiveZScores.append(getZScore(splitLine[i],activeMeans[i-1],activeSDs[i-1]))
            myOutString += 'activeZScores\t'+joiner(myActiveZScores)+'\n'
            myInactiveZScores = []
            for i in range(1,len(myHeader)):
                myInactiveZScores.append(getZScore(splitLine[i],inactiveMeans[i-1],inactiveSDs[i-1]))
            myOutString += 'inactiveZScores\t'+joiner(myInactiveZScores)+'\n'
            myOutString += 'lowerZScore\t'+compareZScores(myActiveZScores,myInactiveZScores)+'\n\n'

    open('mouseOutliers.txt','w').write(myOutString)


def getZScore(val,m,sd):
    if not str(val) == '?':
        return((float(val)-m)/sd)
    else:
        return('-')

def compareZScores(active,inactive):
    myReturn = []
    for i in range(0,len(active)):
        if not str(active[i]) == '-':
            if abs(float(active[i])) < abs(float(inactive[i])):
                myReturn.append('active')
            else:
                myReturn.append('inactive')
        else:
            myReturn.append('-')
    return(joiner(myReturn))




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
    getMouseOutlierData()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

