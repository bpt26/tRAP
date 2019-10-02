#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# parsePredictions.py

import sys
import os
import time
import random
import numpy
import gzip
import math

def gettRNAs(myAlignment):
    mytRNAs = {}
    for line in open(myAlignment):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'Ananc2':
            for k in splitLine:
                if not k == 'N/A':
                    mytRNAs[k] = True
    return(mytRNAs)



def getDistribs(abb, alignedtRNAs):
    isoTypeToCounts = {}
    myKeys = ['AlaAGC','AlaCGC','AlaGGC','AlaTGC','ArgACG','ArgCCG','ArgCCT','ArgGCG','ArgTCG']
    myKeys += ['ArgTCT','AsnATT','AsnGTT','AspATC','AspGTC','CysACA','CysGCA','GlnCTG','GlnTTG']
    myKeys += ['GluCTC','GluTTC','GlyACC','GlyCCC','GlyGCC','GlyTCC','HisATG','HisGTG','IleAAT']
    myKeys += ['IleGAT','IleTAT','LeuAAG','LeuCAA','LeuCAG','LeuTAA','LeuTAG','LysCTT','LysTTT']
    myKeys += ['MetCAT','PheAAA','PheGAA','ProAGG','ProCGG','ProGGG','ProTGG','SeCTCA','SerACT']
    myKeys += ['SerAGA','SerCGA','SerGCT','SerGGA','SerTGA','SupCTA','SupTCA','SupTTA','ThrAGT']
    myKeys += ['ThrGGT','ThrCGT','ThrTGT','TrpCCA','TyrATA','TyrGTA','UndetNNN','ValAAC','ValCAC','ValGAC','ValTAC','iMetCAT']
    for k in myKeys:
        isoTypeToCounts[k] = [0, 0]

    tRNATotal = 0
    tRNAActive = 0
    tRNAInactive = 0
    myOut = 'anticodonCountBySpeciesNewNoSegDups.txt'
    mySpecies = abb
    tRNATotal = 0

    for line in open(abb.upper()+'/'+abb+'tRNAClassificationsNewNoSegDups.txt'):
        splitLine = (line.strip()).split('\t')
        mytRNAName = splitLine[0]
        if abb+'-'+mytRNAName in alignedtRNAs:
            if not abb == 'Hsapi38':
                myIso = mytRNAName.split('-')[-1]
            else:
                myIso = ''.join(mytRNAName.split('-')[1:3])
            myConf = float(splitLine[1])
            tRNATotal += 1
            if myConf > 0:
                (isoTypeToCounts[myIso][0]) += 1
                tRNAActive += 1
            else:
                (isoTypeToCounts[myIso][1]) += 1
                tRNAInactive += 1

    if mySpecies == 'Ananc2':
        myOutString = 'species\t'+joiner(sorted(myKeys))+'\ttotal\tactive\tinactive\n'
    else:
        myOutString = ''
    myOutString += mySpecies
    for k in sorted(myKeys):
        myOutString += '\t'+joinerC(isoTypeToCounts[k])
    myOutString += '\t'+str(tRNATotal)+'\t'+str(tRNAActive)+'\t'+str(tRNAInactive)+'\n'
    print(mySpecies, tRNAActive, tRNAInactive)
    if mySpecies == 'Ananc2':
        open(myOut, 'w').write(myOutString)
    else:
        open(myOut, 'a').write(myOutString)
    sys.stderr.write("Finished "+abb+'\n')

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def joinerC(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ','.join(newList)



def main():

    alignedtRNAs = gettRNAs('oneToOneMapAllSpeciesAugmentMacaqueNoSegDups.txt')
    mySpeciesAbbrevs = ['Hsapi38','Ptrog5','Ggori5','Ppygm3','Mmula8','Ananc2','Mmuri3','Jjacu1','Mochr1','Mmusc10']
    mySpeciesAbbrevs += ['Rnorv6','Hglab2','Cporc3','Clani1','Odegu1','Mmarm2','Ocuni2','Sscro11','Oorca1','Btaur8']
    mySpeciesAbbrevs += ['Chirc1','Ecaba2','Cfami3','Mputo1','Palec1','Efusc1','Eeuro2','Casia1','Dnove3']
    for k in sorted(mySpeciesAbbrevs):
        getDistribs(k, alignedtRNAs)

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit


