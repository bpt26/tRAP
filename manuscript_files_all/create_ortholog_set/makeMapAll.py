#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 10/15/2018
# makePrimateMap.py

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

def makePrimateMap():

    myWorkingDir = '/public/groups/corbettlab/tRNA/classifier/JOEL/'

    abbrevToFull = {}
    abbrevToFull['Homo_sapiens'] = 'Hsapi38'
    abbrevToFull['Hsapi38'] = 'Homo_sapiens'
    abbrevToFull['Pan_troglodytes'] = 'Ptrog5'
    abbrevToFull['Ptrog5'] = 'Pan_troglodytes'
    abbrevToFull['Gorilla_gorilla_gorilla'] = 'Ggori5'
    abbrevToFull['Ggori5'] = 'Gorilla_gorilla_gorilla'
    abbrevToFull['Pongo_abelii'] = 'Ppygm3'
    abbrevToFull['Ppygm3'] = 'Pongo_abelii'
    abbrevToFull['Macaca_mulatta'] = 'Mmula8'
    abbrevToFull['Mmula8'] = 'Macaca_mulatta'
    abbrevToFull['Aotus_nancymaae'] = 'Ananc2'
    abbrevToFull['Ananc2'] = 'Aotus_nancymaae'
    abbrevToFull['Microcebus_murinus'] = 'Mmuri3'
    abbrevToFull['Mmuri3'] = 'Microcebus_murinus'
    abbrevToFull['Jaculus_jaculus'] = 'Jjacu1'
    abbrevToFull['Jjacu1'] = 'Jaculus_jaculus'
    abbrevToFull['Microtus_ochrogaster'] = 'Mochr1'
    abbrevToFull['Mochr1'] = 'Microtus_ochrogaster'
    abbrevToFull['Mus_musculus'] = 'Mmusc10'
    abbrevToFull['Mmusc10'] = 'Mus_musculus'
    abbrevToFull['Rattus_norvegicus'] = 'Rnorv6'
    abbrevToFull['Rnorv6'] = 'Rattus_norvegicus'
    abbrevToFull['Heterocephalus_glaber'] = 'Hglab2'
    abbrevToFull['Hglab2'] = 'Heterocephalus_glaber'
    abbrevToFull['Cavia_porcellus'] = 'Cporc3'
    abbrevToFull['Cporc3'] = 'Cavia_porcellus'
    abbrevToFull['Chinchilla_lanigera'] = 'Clani1'
    abbrevToFull['Clani1'] = 'Chinchilla_lanigera'
    abbrevToFull['Octodon_degus'] = 'Odegu1'
    abbrevToFull['Odegu1'] = 'Octodon_degus'
    abbrevToFull['Marmota_marmota'] = 'Mmarm2'
    abbrevToFull['Mmarm2'] = 'Marmota_marmota'
    abbrevToFull['Oryctolagus_cuniculus'] = 'Ocuni2'
    abbrevToFull['Ocuni2'] = 'Oryctolagus_cuniculus'
    abbrevToFull['Sus_scrofa'] = 'Sscro11'
    abbrevToFull['Sscro11'] = 'Sus_scrofa'
    abbrevToFull['Orcinus_orca'] = 'Oorca1'
    abbrevToFull['Oorca1'] = 'Orcinus_orca'
    abbrevToFull['Bos_taurus'] = 'Btaur8'
    abbrevToFull['Btaur8'] = 'Bos_taurus' 
    abbrevToFull['Capra_hircus'] = 'Chirc1'
    abbrevToFull['Chirc1'] = 'Capra_hircus'
    abbrevToFull['Equus_caballus'] = 'Ecaba2'
    abbrevToFull['Ecaba2'] = 'Equus_caballus'
    abbrevToFull['Canis_lupus_familiaris'] = 'Cfami3'
    abbrevToFull['Cfami3'] = 'Canis_lupus_familiaris'
    abbrevToFull['Mustela_putorius_furo'] = 'Mputo1'
    abbrevToFull['Mputo1'] = 'Mustela_putorius_furo'
    abbrevToFull['Pteropus_alecto'] = 'Palec1'
    abbrevToFull['Palec1'] = 'Pteropus_alecto'
    abbrevToFull['Eptesicus_fuscus'] = 'Efusc1'
    abbrevToFull['Efusc1'] = 'Eptesicus_fuscus'
    abbrevToFull['Erinaceus_europaeus'] = 'Eeuro2'
    abbrevToFull['Eeuro2'] = 'Erinaceus_europaeus'
    abbrevToFull['Chrysochloris_asiatica'] = 'Casia1'
    abbrevToFull['Casia1'] = 'Chrysochloris_asiatica'
    abbrevToFull['Dasypus_novemcinctus'] = 'Dnove3'
    abbrevToFull['Dnove3'] = 'Dasypus_novemcinctus'

    mySpeciesNames = []
    mySpeciesAbbrevs = []
    for k in (abbrevToFull.keys()):
        if '_' in k:
            mySpeciesNames.append(k)
        else:
            mySpeciesAbbrevs.append(k)

    tRNATotRNAToAlignCount = {}

    coordTotRNA = {}
    tRNAToCoords = {}
    halChromToChr = {}

    for species in mySpeciesAbbrevs:
        myAllBed = species+'_tRNAs.bed'
        for line in open(myWorkingDir+myAllBed):
            splitLine = (line.strip()).split('\t')
            tRNAToCoords[species+'-'+splitLine[3]] = []
            tRNATotRNAToAlignCount[species+'-'+splitLine[3]] = {}
            for k in range(int(splitLine[1]), int(splitLine[2])):
                coordTotRNA[species+'-'+str(splitLine[0])+'-'+str(k)] = species+'-'+splitLine[3]
                (tRNAToCoords[species+'-'+splitLine[3]]).append(species+'-'+str(splitLine[0])+'-'+str(k))
            if 'chr' in splitLine[3]:
                print(species+'-'+splitLine[3])
                print(tRNAToCoords[species+'-'+splitLine[3]])
    print(len(tRNAToCoords))
    humanCount = 0
    for k in tRNAToCoords:
        if k.startswith('Hsapi38-'):
            humanCount += 1
    print(humanCount)


    for species in sorted(mySpeciesAbbrevs):
        myMAF = species+'_tRNAs.maf'
        for line in open(myWorkingDir+myMAF):
            if len(line.strip()) > 0:
                splitLine = (line.strip()).split('\t')
                if splitLine[0] == 'a':
                    prevLine = 'a'
                if splitLine[0] == 's':
                    mySpecies = (splitLine[1]).split('.')[0]
                    myChrom = (splitLine[1]).split('.')[1]
                    if len((splitLine[1]).split('.')) > 2:
                        myChrom = '.'.join((splitLine[1]).split('.')[1:])
                    if mySpecies == abbrevToFull[species] and prevLine == 'a':
                        myCurrenttRNAs = []
                        for k in range(int(splitLine[2]), int(splitLine[2])+int(splitLine[3])):
                            if abbrevToFull[mySpecies]+'-'+myChrom+'-'+str(k) in coordTotRNA:
                                myCurrenttRNA = coordTotRNA[abbrevToFull[mySpecies]+'-'+myChrom+'-'+str(k)]
                                if not myCurrenttRNA in myCurrenttRNAs:
                                    myCurrenttRNAs.append(myCurrenttRNA)
                        if len(myCurrenttRNAs) > 0:
                            myCurrenttRNA = myCurrenttRNAs[0]
                        else:
                            myCurrenttRNA = ''
                        myCurrenttRNAs = []
                        prevLine = ''
                    elif mySpecies in abbrevToFull:
                        if not myCurrenttRNA == '':
                            for k in range(int(splitLine[2]), int(splitLine[2])+int(splitLine[3])):
                                if abbrevToFull[mySpecies]+'-'+myChrom+'-'+str(k) in coordTotRNA:
                                    mytRNA = coordTotRNA[abbrevToFull[mySpecies]+'-'+myChrom+'-'+str(k)]
                                    if not mytRNA in tRNATotRNAToAlignCount[myCurrenttRNA]:
                                        (tRNATotRNAToAlignCount[myCurrenttRNA])[mytRNA] = 0.0
                                    (tRNATotRNAToAlignCount[myCurrenttRNA])[mytRNA] += 1
                        prevLine = ''
        sys.stderr.write('Finished reading alignments for '+species+'\n')

    myOut = 'cactustRNAOrthologMap.txt'
    myString = ''
    for tRNA in sorted(tRNATotRNAToAlignCount.keys()):
        myString += tRNA
        for align in tRNATotRNAToAlignCount[tRNA]:
            myString += '\t'+align+'\t'+str((tRNATotRNAToAlignCount[tRNA])[align])
        myString += '\n'
    open(myOut, 'w').write(myString)

    myCount1 = 0
    myCount2 = 0
    for tRNA in sorted(tRNATotRNAToAlignCount.keys()):
        if tRNA.startswith('Hsapi38'):
            myCount2 += 1
            if 'chr' in tRNA:
                print(tRNA)
    #print(myCount1)
    print(myCount2)


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
    makePrimateMap()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit










