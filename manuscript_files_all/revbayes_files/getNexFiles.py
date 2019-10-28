#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 11/29/2018
# getPrimateSet.py

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

def getNexFiles():
    abbrevToFull = {}
    abbrevToFull['Hsapi38'] = 'Homo_sapiens'
    abbrevToFull['Ptrog5'] = 'Pan_troglodytes'
    abbrevToFull['Ggori5'] = 'Gorilla_gorilla_gorilla'
    abbrevToFull['Ppygm3'] = 'Pongo_abelii'
    abbrevToFull['Mmula8'] = 'Macaca_mulatta'
    abbrevToFull['Ananc2'] = 'Aotus_nancymaae'
    abbrevToFull['Mmuri3'] = 'Microcebus_murinus'
    abbrevToFull['Jjacu1'] = 'Jaculus_jaculus'
    abbrevToFull['Mochr1'] = 'Microtus_ochrogaster'
    abbrevToFull['Mmusc10'] = 'Mus_musculus'
    abbrevToFull['Rnorv6'] = 'Rattus_norvegicus'
    abbrevToFull['Hglab2'] = 'Heterocephalus_glaber'
    abbrevToFull['Cporc3'] = 'Cavia_porcellus'
    abbrevToFull['Clani1'] = 'Chinchilla_lanigera'
    abbrevToFull['Odegu1'] = 'Octodon_degus'
    abbrevToFull['Mmarm2'] = 'Marmota_marmota'
    abbrevToFull['Ocuni2'] = 'Oryctolagus_cuniculus'
    abbrevToFull['Sscro11'] = 'Sus_scrofa'
    abbrevToFull['Oorca1'] = 'Orcinus_orca'
    abbrevToFull['Btaur8'] = 'Bos_taurus' 
    abbrevToFull['Chirc1'] = 'Capra_hircus'
    abbrevToFull['Ecaba2'] = 'Equus_caballus'
    abbrevToFull['Cfami3'] = 'Canis_lupus_familiaris'
    abbrevToFull['Mputo1'] = 'Mustela_putorius_furo'
    abbrevToFull['Palec1'] = 'Pteropus_alecto'
    abbrevToFull['Efusc1'] = 'Eptesicus_fuscus'
    abbrevToFull['Eeuro2'] = 'Erinaceus_europaeus'
    abbrevToFull['Casia1'] = 'Chrysochloris_asiatica'
    abbrevToFull['Dnove3'] = 'Dasypus_novemcinctus'

    mySpeciesAbbrevs = ['Hsapi38','Ptrog5','Ggori5','Ppygm3','Mmula8','Ananc2','Mmuri3','Jjacu1','Mochr1','Mmusc10']
    mySpeciesAbbrevs += ['Rnorv6','Hglab2','Cporc3','Clani1','Odegu1','Mmarm2','Ocuni2','Sscro11','Oorca1','Btaur8']
    mySpeciesAbbrevs += ['Chirc1','Ecaba2','Cfami3','Mputo1','Palec1','Efusc1','Eeuro2','Casia1','Dnove3']
    myAbbSort = sorted(mySpeciesAbbrevs)

    tRNAToActivity = {}
    abbToNexString = {}

    for k in myAbbSort:
        abbToNexString[k] = ''
        if not k == 'Hsapi38':
            for line in open(k.upper()+'/'+k+'tRNAClassificationsNewFixedWithSegDups.txt'):
                splitLine = (line.strip()).split('\t')
                if splitLine[2] == 'active':
                    tRNAToActivity[k+'-'+splitLine[0]] = '2'
                elif splitLine[2] == 'inactive':
                    tRNAToActivity[k+'-'+splitLine[0]] = '1'
    for line in open('HSAPI38/allHumantRNAPredictions.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[2] == 'active':
            tRNAToActivity['Hsapi38-'+splitLine[0]] = '2'
        else:
            tRNAToActivity['Hsapi38-'+splitLine[0]] = '1'

    primateToString = {}
    clade2ToString = {}
    clade3ToString = {}

    for line in open('FINAL_ORTHOLOG_SET.txt'):
        splitLine = (line.strip()).split('\t')

        for i in range(0,len(splitLine)):
            if (splitLine[13] in tRNAToActivity or splitLine[13] == 'N/A') and (splitLine[18] in tRNAToActivity or splitLine[18] == 'N/A'):
                if splitLine[i] in tRNAToActivity:
                    abbToNexString[myAbbSort[i]] += tRNAToActivity[splitLine[i]]
                elif splitLine[i] == 'N/A':
                    abbToNexString[myAbbSort[i]] += '0'
                else:
                    print(splitLine[i])

        if (splitLine[13] in tRNAToActivity or splitLine[13] == 'N/A'):
            if checkMyLine(splitLine,[0,11,13,16,17,25,26],tRNAToActivity):
                for i in [0,11,13,16,17,25,26]:
                    if not myAbbSort[i] in primateToString:
                        primateToString[myAbbSort[i]] = ''
                    if splitLine[i] in tRNAToActivity:
                        primateToString[myAbbSort[i]] += tRNAToActivity[splitLine[i]]
                    elif splitLine[i] == 'N/A':
                        primateToString[myAbbSort[i]] += '0'

        if (splitLine[18] in tRNAToActivity or splitLine[18] == 'N/A'):
            if checkMyLine(splitLine,[5,6,12,14,15,18,19,21,22,27],tRNAToActivity):
                for i in [5,6,12,14,15,18,19,21,22,27]:
                    if not myAbbSort[i] in clade2ToString:
                        clade2ToString[myAbbSort[i]] = ''
                    if splitLine[i] in tRNAToActivity:
                        clade2ToString[myAbbSort[i]] += tRNAToActivity[splitLine[i]]
                    elif splitLine[i] == 'N/A':
                        clade2ToString[myAbbSort[i]] += '0'

        if checkMyLine(splitLine,[1,3,4,8,9,10,20,23,24,28],tRNAToActivity):
            for i in [1,3,4,8,9,10,20,23,24,28]:
                if not myAbbSort[i] in clade3ToString:
                    clade3ToString[myAbbSort[i]] = ''
                if splitLine[i] in tRNAToActivity:
                    clade3ToString[myAbbSort[i]] += tRNAToActivity[splitLine[i]]
                elif splitLine[i] == 'N/A':
                    clade3ToString[myAbbSort[i]] += '0'

    myAllSpeciesOutString = '#NEXUS\nBEGIN DATA;\n\tDIMENSIONS  NCHAR='+str(len(abbToNexString['Ananc2']))+' NTAX=29;\n\tFORMAT '
    myAllSpeciesOutString += 'DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="012";\n\n\tMATRIX\n'
    primateOutString = '#NEXUS\nBEGIN DATA;\n\tDIMENSIONS  NCHAR='+str(len(primateToString['Ptrog5']))+' NTAX=7;\n\tFORMAT '
    primateOutString += 'DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="012";\n\n\tMATRIX\n'
    clade2OutString = '#NEXUS\nBEGIN DATA;\n\tDIMENSIONS  NCHAR='+str(len(clade2ToString['Mmusc10']))+' NTAX=10;\n\tFORMAT '
    clade2OutString += 'DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="012";\n\n\tMATRIX\n'
    clade3OutString = '#NEXUS\nBEGIN DATA;\n\tDIMENSIONS  NCHAR='+str(len(clade3ToString['Btaur8']))+' NTAX=10;\n\tFORMAT '
    clade3OutString += 'DATATYPE=STANDARD GAP=- MISSING=? SYMBOLS="012";\n\n\tMATRIX\n'
    for k in myAbbSort:
        myAllSpeciesOutString += '\t'+abbrevToFull[k]+'\t'+abbToNexString[k]+'\n'
        if k in primateToString:
            primateOutString += '\t'+abbrevToFull[k]+'\t'+primateToString[k]+'\n'
        if k in clade2ToString:
            clade2OutString += '\t'+abbrevToFull[k]+'\t'+clade2ToString[k]+'\n'
        if k in clade3ToString:
            clade3OutString += '\t'+abbrevToFull[k]+'\t'+clade3ToString[k]+'\n'

    myAllSpeciesOutString += '\n;\n\nEND;\n'
    primateOutString += '\n;\n\nEND;\n'
    clade2OutString += '\n;\n\nEND;\n'
    clade3OutString += '\n;\n\nEND;\n'

    open('allSpeciesMatrix.nex', 'w').write(myAllSpeciesOutString)
    open('primateMatrix.nex', 'w').write(primateOutString)
    open('clade2Matrix.nex', 'w').write(clade2OutString)
    open('clade3Matrix.nex', 'w').write(clade3OutString)


def checkMyLine(myLine, myIndices, myDict):
    for k in myIndices:
        if myLine[k] in myDict:
            return(True)
    return(False)


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
    getNexFiles()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit