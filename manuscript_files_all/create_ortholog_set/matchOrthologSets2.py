#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# matchOrthologSets.py

import sys
import os.path
import time
import random
import numpy
import gzip
import math

##########################
##### MAIN FUNCTIONS #####
##########################

def matchOrthologSetsCow():
    bt8CoordsTotRNA = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Btaur8/bosTau8-tRNAsnew.bed'):
        splitLine = (line.strip()).split('\t')
        for k in range(int(splitLine[1]),int(splitLine[2])):
            bt8CoordsTotRNA[splitLine[0]+'-'+str(k)] = splitLine[3]

    halChromToRealChrom = {}
    realChromToHalChrom = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Btaur8/Btaur8_assembly_report.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 6 and splitLine[1] == 'assembled-molecule':
            halChromToRealChrom[splitLine[6]] = (splitLine[0]).lower()
            realChromToHalChrom[(splitLine[0]).lower()] = splitLine[6]

    bt8NameToHalName = {}
    halNameToBt8Name = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Btaur8/Btaur8_tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        if halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25) in bt8CoordsTotRNA:
            halName = splitLine[3]
            bt8Name = bt8CoordsTotRNA[halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25)]
            bt8NameToHalName[bt8Name] = halName
            halNameToBt8Name[halName] = bt8Name

    bt8ToHuman = {}
    humanToBt8 = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/andrewOrthologs.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[4] in ['None', 'hg19-bosTau8-pairs.txt']:
            bt8ToHuman[splitLine[4]] = splitLine[0]
            humanToBt8[splitLine[0]] = splitLine[4]
    humanToMaps = {}
    cowToMaps = {}
    myOut = 'oneToOneMapAllSpeciesAugmentCow.txt'
    myOutString = ''
    myOutPreCombine = ''
    for line in open('oneToOneMapAllSpeciesAugment.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'Ananc2':
            myHeader = splitLine
            myOutString += joiner(splitLine)+'\n'
            myOutPreCombine += joiner(splitLine)+'\n'
        else:
            if not splitLine[myHeader.index('Hsapi38')] == 'N/A':
                humanToMaps[splitLine[myHeader.index('Hsapi38')].split('Hsapi38-')[-1]] = splitLine
            if not splitLine[myHeader.index('Btaur8')] == 'N/A':
                cowToMaps[splitLine[myHeader.index('Btaur8')].split('Btaur8-')[-1]] = splitLine

    count = 0
    dontPrint = [myHeader]
    for tRNA in humanToBt8:
        if tRNA in humanToMaps:
            print('humantRNA: '+tRNA+'\t'+joiner(humanToMaps[tRNA]))
            humanSpecies = getSpecies(humanToMaps[tRNA])
            if humanToBt8[tRNA] in bt8NameToHalName:
                myCowtRNA = bt8NameToHalName[humanToBt8[tRNA]]
                if myCowtRNA in cowToMaps:
                    cowSpecies = getSpecies(cowToMaps[myCowtRNA])
                    print('cowtRNA: '+myCowtRNA+'\t'+joiner(cowToMaps[myCowtRNA]))
                if humanToMaps[tRNA] == cowToMaps[myCowtRNA]:
                    count += 1
                if getOverlap(cowSpecies,humanSpecies):
                    myOutString += combineLists(cowToMaps[myCowtRNA],humanToMaps[tRNA])
                    myOutPreCombine += joiner(humanToMaps[tRNA])+'\n'+joiner(cowToMaps[myCowtRNA])+'\n'
                    dontPrint.append(cowToMaps[myCowtRNA])
                    dontPrint.append(humanToMaps[tRNA])

    open('onlyCombinedLinesCow.txt', 'w').write(myOutString)
    open('beforeCombiningLinesCow.txt', 'w').write(myOutPreCombine)
    for line in open('oneToOneMapAllSpeciesAugment.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine in dontPrint:
            myOutString += joiner(splitLine)+'\n'
    open(myOut, 'w').write(myOutString)



    print(count)
    print(len(dontPrint))




def joiner(entry):
    """
    Helper function to print lists in 
    tab separated format.
    """
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def getSpecies(myList):
    """
    Get list of species not N/A in ortholog set.
    """
    toReturn = []
    for k in myList:
        if '-' in k:
            toReturn.append(k.split('-')[0])
    return(toReturn)

def getOverlap(list1, list2):
    """
    Return false if any species shared by two lists, else true.
    """
    for k in list1:
        if k in list2:
            return(False)
    return(True)

def combineLists(list1,list2):
    """
    Return combined list of non-lappping input lists.
    """
    toReturn = []
    for i in range(0,len(list1)):
        if list1[i] != 'N/A':
            toReturn.append(list1[i])
        else:
            toReturn.append(list2[i])
        i += 1
    return(joiner(toReturn)+'\n')

def main():
    matchOrthologSetsCow()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit