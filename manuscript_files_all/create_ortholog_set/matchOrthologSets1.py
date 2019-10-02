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

def matchOrthologSets():
    mm10CoordsTotRNA = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmusc10/mm10-tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        for k in range(int(splitLine[1]),int(splitLine[2])):
            mm10CoordsTotRNA[splitLine[0]+'-'+str(k)] = splitLine[3]

    halChromToRealChrom = {}
    realChromToHalChrom = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmusc10/Mmusc10_assembly_report.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 6 and splitLine[1] == 'assembled-molecule':
            halChromToRealChrom[splitLine[6]] = splitLine[9]
            realChromToHalChrom[splitLine[9]] = splitLine[6]

    mm10NameToHalName = {}
    halNameToMm10Name = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/Mmusc10/Mmusc10_tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        if halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25) in mm10CoordsTotRNA:
            halName = splitLine[3]
            mm10Name = mm10CoordsTotRNA[halChromToRealChrom[splitLine[0]]+'-'+str(int(splitLine[1])+25)]
            mm10NameToHalName[mm10Name] = halName
            halNameToMm10Name[halName] = mm10Name

    mm10ToHuman = {}
    humanToMm10 = {}
    bt8ToHuman = {}
    humanToBt8 = {}
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/andrewOrthologs.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[1] in ['None', 'hg19-mm10-pairs.txt']:
            mm10ToHuman[splitLine[1]] = splitLine[0]
            humanToMm10[splitLine[0]] = splitLine[1]

    mouseToMaps = {}
    humanToMaps = {}
    cowToMaps = {}
    myOut = 'oneToOneMapAllSpeciesAugment.txt'
    myOutString = ''
    myOutPreCombine = ''
    myOutSame = ''
    myOutDisagree = ''
    for line in open('oneToOneMapAllSpecies.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'Ananc2':
            myHeader = splitLine
            myOutString += joiner(splitLine)+'\n'
            myOutPreCombine += joiner(splitLine)+'\n'
        else:
            if not splitLine[myHeader.index('Mmusc10')] == 'N/A':
                mouseToMaps[splitLine[myHeader.index('Mmusc10')].split('Mmusc10-')[-1]] = splitLine
            if not splitLine[myHeader.index('Hsapi38')] == 'N/A':
                humanToMaps[splitLine[myHeader.index('Hsapi38')].split('Hsapi38-')[-1]] = splitLine

    for k in humanToMaps:
        print(k, humanToMaps[k])
    for k in mouseToMaps:
        print(k, mouseToMaps[k])
    count = 0
    same = 0
    disagree = 0
    exact = 0
    dontPrint = [myHeader]
    for tRNA in humanToMm10:
        if tRNA in humanToMaps:
            humanSpecies = getSpecies(humanToMaps[tRNA])
            myMousetRNA = mm10NameToHalName[humanToMm10[tRNA]]
            print(myMousetRNA)
            if myMousetRNA in mouseToMaps:
                mouseSpecies = getSpecies(mouseToMaps[myMousetRNA])
                if humanToMaps[tRNA] == mouseToMaps[myMousetRNA]:
                    exact += 1
                

            if getOverlap(mouseSpecies,humanSpecies):
                myOutString += combineLists(mouseToMaps[myMousetRNA],humanToMaps[tRNA])
                myOutPreCombine += joiner(humanToMaps[tRNA])+'\n'+joiner(mouseToMaps[myMousetRNA])+'\n'
                dontPrint.append(mouseToMaps[myMousetRNA])
                dontPrint.append(humanToMaps[tRNA])
                count += 1
            elif sorted(mouseSpecies) == sorted(humanSpecies):
                myOutSame += joiner(humanToMaps[tRNA])+'\n'+joiner(mouseToMaps[myMousetRNA])+'\n'
                same += 1
            else:
                myOutDisagree += joiner(humanToMaps[tRNA])+'\n'+joiner(mouseToMaps[myMousetRNA])+'\n'
                disagree += 1


    open('onlyCombinedLines.txt', 'w').write(myOutString)
    open('beforeCombiningLines.txt', 'w').write(myOutPreCombine)
    open('outSame.txt', 'w').write(myOutSame)
    open('outDisagree.txt', 'w').write(myOutDisagree)
    for line in open('oneToOneMapAllSpecies.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine in dontPrint:
            myOutString += joiner(splitLine)+'\n'
    open(myOut, 'w').write(myOutString)

    print(exact)
    print(count)
    print(same)
    print(disagree)


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
    matchOrthologSets()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit
