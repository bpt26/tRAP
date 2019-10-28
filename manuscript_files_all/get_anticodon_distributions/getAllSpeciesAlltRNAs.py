#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 2/1/2018
# compareDatabases.py

import sys
import os
import datetime
import random
import numpy
import gzip
import math


##########################
##### MAIN FUNCTIONS #####
##########################

# MANUALLY CHANGE KV888377.1.tRNA1-LysTTT to tRNA-Lys-TTT-1-7

def combineAllTables():
    mySpecies = ['ANANC2','BTAUR8','CASIA1','CFAMI3','CHIRC1','CLANI1','CPORC3','DNOVE3','ECABA2','EEURO2','EFUSC1','GGORI5']
    mySpecies += ['HGLAB2','HSAPI38','JJACU1','MMARM2','MMULA8','MMURI3','MMUSC10','MOCHR1','MPUTO1','OCUNI2','ODEGU1','OORCA1']
    mySpecies += ['PALEC1','PPYGM3','PTROG5','RNORV6','SSCRO11']

    speciesToFull = {}
    speciesToFull['ANANC2'] = 'Aotus_nancymaae'
    speciesToFull['BTAUR8'] = 'Bos_taurus'
    speciesToFull['CASIA1'] = 'Chrysochloris_asiatica' 
    speciesToFull['CFAMI3'] = 'Canis_lupus_familiaris'
    speciesToFull['CHIRC1'] = 'Capra_hircus'
    speciesToFull['CLANI1'] = 'Chinchilla_lanigera'
    speciesToFull['CPORC3'] = 'Cavia_porcellus'
    speciesToFull['DNOVE3'] = 'Dasypus_novemcinctus'
    speciesToFull['ECABA2'] = 'Equus_caballus'
    speciesToFull['EEURO2'] = 'Erinaceus_europaeus'
    speciesToFull['EFUSC1'] = 'Eptesicus_fuscus'
    speciesToFull['GGORI5'] = 'Gorilla_gorilla'
    speciesToFull['HGLAB2'] = 'Heterocephalus_glaber'
    speciesToFull['HSAPI38'] = 'Homo_sapiens'
    speciesToFull['JJACU1'] = 'Jaculus_jaculus'
    speciesToFull['MMARM2'] = 'Marmota_marmota'
    speciesToFull['MMULA8'] = 'Macaca_mulatta'
    speciesToFull['MMURI3'] = 'Microcebus_murinus'
    speciesToFull['MMUSC10'] = 'Mus_musculus'
    speciesToFull['MOCHR1'] = 'Microtus_ochrogaster'
    speciesToFull['MPUTO1'] = 'Mustela_putorius'
    speciesToFull['OCUNI2'] = 'Oryctolagus_cuniculus'
    speciesToFull['ODEGU1'] = 'Octodon_degus'
    speciesToFull['OORCA1'] = 'Orcinus_orca'
    speciesToFull['PALEC1'] = 'Pteropus_alecto'
    speciesToFull['PPYGM3'] = 'Pongo_abelii'
    speciesToFull['PTROG5'] = 'Pan_troglodytes'
    speciesToFull['RNORV6'] = 'Rattus_norvegicus'
    speciesToFull['SSCRO11'] = 'Sus_scrofa'

    myHeader = 'species\ttRNA\tchr\tstart\tend\tstrand\t'
    #myHeader += 'tRNAPhyloPAvg\t5PhyloPAvg\tCpGDensity\tObsExpUp\tGenBit\ttRNA10kb\tProt75kb\tTTTT\tCodon\tMFE\tprediction\tprobability\tinTestSet\n'
    myHeader += 'Total Number of tRNA Genes With Identical Anticodon\tMinimum Free Energy of Canonical tRNA Secondary Structure\ttRNAscan-SE General Bit Score\t'
    myHeader += 'Average PhyloP Score in tRNA Gene Sequence\tDistance To Nearest TTTT Transcription Termination Sequence\tCpG Density Across tRNA Locus\t'
    myHeader += "Observed/Expected CpG Islands Score Upstream of tRNA Gene\tAverage PhyloP Score in 5' Flanking Region\t"
    myHeader += 'tRNA Genes Within 10 Kilobases\tExons Within 75 Kilobases\tprediction\tprobability\tinTestSet\n'
    myOutString = myHeader
    for species in mySpecies:
        tRNAToChrStartEndStrand = {}
        for line in open(species+'/tRNAHiConf.bed'):
            splitLine = (line.strip()).split('\t')
            tRNAToChrStartEndStrand[splitLine[3]] = [splitLine[0],splitLine[1],splitLine[2],splitLine[5]]
        tRNAToData = {} 
        for line in open(species+'/alltRNAData.tsv'):
            s = (line.strip()).split('\t')
            if s[0] != 'tRNA':
                tRNAToData[s[0]] = [s[9],s[10],s[5],s[1],s[8],s[3],s[4],s[2],s[6],s[7]]
        tRNAToTestSet = {}
        tRNAToPred = {}
        mytRNAs = []
        for line in open(species+'/'+partLower(species)+'tRNAClassificationsNewFixedNoSegDups.txt'):
            splitLine = (line.strip()).split('\t')
            tRNAToTestSet[splitLine[0]] = True
        for line in open(species+'/'+partLower(species)+'tRNAClassificationsNewWithSegDups.txt'):
            splitLine = (line.strip()).split('\t')
            tRNAToPred[splitLine[0]] = [splitLine[2],splitLine[1]]
            mytRNAs.append(splitLine[0])
        print(species)
        for t in sorted(mytRNAs):
            myOutString += speciesToFull[species]+'\t'+t+'\t'+joiner(tRNAToChrStartEndStrand[t])+'\t'+joiner(tRNAToData[t])+'\t'+joiner(tRNAToPred[t])+'\t'
            if t in tRNAToTestSet:
                myOutString += 'True\n'
            else:
                myOutString += 'False\n'
    open('all_species_all_trnas_fixed.txt','w').write(myOutString)


def joiner(entry):
    """
    Helper function to print lists in 
    tab separated format.
    """
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def partLower(myString):
    return(myString[0]+(myString[1:].lower()))




def main():
    combineAllTables()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

