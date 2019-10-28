#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# checkCoretRNAActivityLevels.py

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

def checkCoretRNAActivityLevels():
    tRNAToActivity = {}
    abbs = ['Ananc2','Btaur8','Casia1','Cfami3','Chirc1','Clani1','Cporc3','Dnove3','Ecaba2','Eeuro2','Efusc1']
    abbs += ['Ggori5','Hglab2','Hsapi38','Jjacu1','Mmarm2','Mmula8','Mmuri3','Mmusc10','Mochr1','Mputo1','Ocuni2']
    abbs += ['Odegu1','Oorca1','Palec1','Ppygm3','Ptrog5','Rnorv6','Sscro11']
    for k in abbs:
        for line in open(k.upper()+'/'+k+'tRNAClassificationsNewFixedNoSegDups.txt'):
            splitLine = (line.strip()).split()
            tRNAToActivity[k+'-'+splitLine[0]] = splitLine[2]

    activityLists = []
    activeCount = 0
    inactiveCount = 0
    for line in open('out37.txt'):
        splitLine = (line.strip()).split('\t')
        tempList = []
        if not splitLine[0] == 'Hsapi38':
            for k in splitLine:
                if k in tRNAToActivity:
                    tempList.append(tRNAToActivity[k])
        activityLists.append(tempList)
        activeCount += tempList.count('active')
        inactiveCount += tempList.count('inactive')

    diffCount = 0
    majActive = 0
    majInactive = 0
    for k in activityLists:
        print(k)
        if 'inactive' in k and 'active' in k:
            diffCount += 1
        if k.count('active') > k.count('inactive'):
            majActive += 1
        elif k.count('inactive') > k.count('active'):
            majInactive += 1
    print(diffCount)
    print(majActive)
    print(majInactive)
    print(float(activeCount)/float(activeCount+inactiveCount)*100.0)

def main():
    checkCoretRNAActivityLevels()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit




