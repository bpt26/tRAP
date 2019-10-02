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

def matchOrthologSetsMacaque():
    mm8CoordsTotRNA = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmula8/rheMac8-tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        for k in range(int(splitLine[1]),int(splitLine[2])):
            mm8CoordsTotRNA[splitLine[0]+'-'+str(k)] = splitLine[3]

    halChromToRealChrom = {}
    realChromToHalChrom = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmula8/Mmula8_assembly_report.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 9 and splitLine[9].startswith('chr'):
            halChromToRealChrom[splitLine[6]] = (splitLine[9])
            realChromToHalChrom[(splitLine[9])] = splitLine[6]

    mm8NameToHalName = {}
    halNameToMm8Name = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmula8/Mmula8_tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        if halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25) in mm8CoordsTotRNA:
            halName = splitLine[3]
            mm8Name = mm8CoordsTotRNA[halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25)]
            mm8NameToHalName[mm8Name] = halName
            halNameToMm8Name[halName] = mm8Name

    mm8ToHuman = {}
    humanToMm8 = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmula8/andrewRheMac3.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'intervalname':
            mm8ToHuman[splitLine[2]] = splitLine[1]
            humanToMm8[splitLine[1]] = splitLine[2]
    humanToMaps = {}
    macaqueToMaps = {}
    myOut = 'oneToOneMapAllSpeciesAugmentMacaque.txt'
    myOutString = ''
    for line in open('oneToOneMapAllSpeciesAugmentDog.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'Ananc2':
            myHeader = splitLine
            myOutString += joiner(splitLine)+'\n'
        else:
            if not splitLine[myHeader.index('Hsapi38')] == 'N/A':
                humanToMaps[splitLine[myHeader.index('Hsapi38')].split('Hsapi38-')[-1]] = splitLine
            if not splitLine[myHeader.index('Mmula8')] == 'N/A':
                macaqueToMaps[splitLine[myHeader.index('Mmula8')].split('Mmula8-')[-1]] = splitLine

    count = 0
    dontPrint = []
    for tRNA in humanToMm8:
        if tRNA in humanToMaps:
            print('humantRNA: '+tRNA+'\t'+joiner(humanToMaps[tRNA]))
            humanSpecies = getSpecies(humanToMaps[tRNA])
            if humanToMm8[tRNA] in mm8NameToHalName:
                myMacaquetRNA = mm8NameToHalName[humanToMm8[tRNA]]
                if myMacaquetRNA in macaqueToMaps:
                    macaqueSpecies = getSpecies(macaqueToMaps[myMacaquetRNA])
                    print('macaquetRNA: '+myMacaquetRNA+'\t'+joiner(macaqueToMaps[myMacaquetRNA]))
                if humanToMaps[tRNA] == macaqueToMaps[myMacaquetRNA]:
                    count += 1
                if getOverlap(macaqueSpecies,humanSpecies):
                    myOutString += combineLists(macaqueToMaps[myMacaquetRNA],humanToMaps[tRNA])
                    dontPrint.append(macaqueToMaps[myMacaquetRNA])
                    dontPrint.append(humanToMaps[tRNA])

    for line in open('oneToOneMapAllSpeciesAugmentDog.txt'):
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
    matchOrthologSetsMacaque()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit