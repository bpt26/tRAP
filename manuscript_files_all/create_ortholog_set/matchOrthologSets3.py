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
    cf3CoordsTotRNA = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Cfami3/canFam3-tRNAsnew.bed'):
        splitLine = (line.strip()).split('\t')
        for k in range(int(splitLine[1]),int(splitLine[2])):
            cf3CoordsTotRNA[splitLine[0]+'-'+str(k)] = splitLine[3]

    halChromToRealChrom = {}
    realChromToHalChrom = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Cfami3/Cfami3_assembly_report.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 6 and splitLine[0].startswith('chr'):
            halChromToRealChrom[splitLine[6]] = (splitLine[9]).lower()
            realChromToHalChrom[(splitLine[9]).lower()] = splitLine[6]

    cf3NameToHalName = {}
    halNameToCf3Name = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Cfami3/Cfami3_tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        if halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25) in cf3CoordsTotRNA:
            halName = splitLine[3]
            cf3Name = cf3CoordsTotRNA[halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25)]
            cf3NameToHalName[cf3Name] = halName
            halNameToCf3Name[halName] = cf3Name

    cf3ToHuman = {}
    humanToCf3 = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/andrewOrthologs.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[4] in ['None', 'hg19-canFam3-pairs.txt']:
            cf3ToHuman[splitLine[4]] = splitLine[0]
            humanToCf3[splitLine[0]] = splitLine[4]
    humanToMaps = {}
    cowToMaps = {}
    myOut = 'oneToOneMapAllSpeciesAugmentDog.txt'
    myOutString = ''
    myOutPreCombine = ''
    for line in open('oneToOneMapAllSpeciesAugmentCow.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'Ananc2':
            myHeader = splitLine
            myOutString += joiner(splitLine)+'\n'
            myOutPreCombine += joiner(splitLine)+'\n'
        else:
            if not splitLine[myHeader.index('Hsapi38')] == 'N/A':
                humanToMaps[splitLine[myHeader.index('Hsapi38')].split('Hsapi38-')[-1]] = splitLine
            if not splitLine[myHeader.index('Cfami3')] == 'N/A':
                cowToMaps[splitLine[myHeader.index('Cfami3')].split('Cfami3-')[-1]] = splitLine

    count = 0
    dontPrint = [myHeader]
    for tRNA in humanToCf3:
        if tRNA in humanToMaps:
            print('humantRNA: '+tRNA+'\t'+joiner(humanToMaps[tRNA]))
            humanSpecies = getSpecies(humanToMaps[tRNA])
            if humanToCf3[tRNA] in cf3NameToHalName:
                myCowtRNA = cf3NameToHalName[humanToCf3[tRNA]]
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

    open('onlyCombinedLinesDog.txt', 'w').write(myOutString)
    open('beforeCombiningLinesDog.txt', 'w').write(myOutPreCombine)
    for line in open('oneToOneMapAllSpeciesAugmentCow.txt'):
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