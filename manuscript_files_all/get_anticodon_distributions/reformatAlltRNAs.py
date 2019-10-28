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

    fullToAbb = {}
    fullToAbb['Aotus_nancymaae'] = 'ANANC2'
    fullToAbb['Bos_taurus'] = 'BTAUR8'
    fullToAbb['Chrysochloris_asiatica'] = 'CASIA1' 
    fullToAbb['Canis_lupus_familiaris'] = 'CFAMI3'
    fullToAbb['Capra_hircus'] = 'CHIRC1'
    fullToAbb['Chinchilla_lanigera'] = 'CLANI1'
    fullToAbb['Cavia_porcellus'] = 'CPORC3'
    fullToAbb['Dasypus_novemcinctus'] = 'DNOVE3'
    fullToAbb['Equus_caballus'] = 'ECABA2'
    fullToAbb['Erinaceus_europaeus'] = 'EEURO2'
    fullToAbb['Eptesicus_fuscus'] = 'EFUSC1'
    fullToAbb['Gorilla_gorilla'] = 'GGORI5'
    fullToAbb['Heterocephalus_glaber'] = 'HGLAB2'
    fullToAbb['Homo_sapiens'] = 'HSAPI38'
    fullToAbb['Jaculus_jaculus'] = 'JJACU1'
    fullToAbb['Marmota_marmota'] = 'MMARM2'
    fullToAbb['Macaca_mulatta'] = 'MMULA8'
    fullToAbb['Microcebus_murinus'] = 'MMURI3'
    fullToAbb['Mus_musculus'] = 'MMUSC10'
    fullToAbb['Microtus_ochrogaster'] = 'MOCHR1'
    fullToAbb['Mustela_putorius'] = 'MPUTO1'
    fullToAbb['Oryctolagus_cuniculus'] = 'OCUNI2'
    fullToAbb['Octodon_degus'] = 'ODEGU1'
    fullToAbb['Orcinus_orca'] = 'OORCA1'
    fullToAbb['Pteropus_alecto'] = 'PALEC1'
    fullToAbb['Pongo_abelii'] = 'PPYGM3'
    fullToAbb['Pan_troglodytes'] = 'PTROG5'
    fullToAbb['Rattus_norvegicus'] = 'RNORV6'
    fullToAbb['Sus_scrofa'] = 'SSCRO11'

    speciesToBed = {}
    speciesToBed['Btaur8'] = 'Btaur8/bosTau8-tRNAs.bed'
    speciesToBed['Cfami3'] = 'Cfami3/canFam3-tRNAs.bed'
    speciesToBed['Cporc3'] = 'Cporc3/cavPor3-tRNAs.bed'
    speciesToBed['Dnove3'] = 'Dnove3/dasNov3-tRNAs.bed'
    speciesToBed['Ecaba2'] = 'Ecaba2/equCab2-tRNAs.bed'
    speciesToBed['Eeuro2'] = 'Eeuro2/eriEur2-tRNAs.bed'
    speciesToBed['Mmula8'] = 'Mmula8/rheMac8-tRNAs.bed'
    speciesToBed['Mmuri3'] = 'Mmuri3/micMur3-tRNAs.bed'
    speciesToBed['Mmusc10'] = 'Mmusc10/mm10-tRNAs.bed'
    speciesToBed['Mputo1'] = 'Mputo1/musFur1-tRNAs.bed'
    speciesToBed['Ocuni2'] = 'Ocuni2/oryCun2-tRNAs.bed'
    speciesToBed['Rnorv6'] = 'Rnorv6/rn6-tRNAs.bed'
    speciesToBed['Sscro11'] = 'Sscro11/susScr11-tRNAs.bed'

    toChange = ['Btaur8','Cfami3','Cporc3','Dnove3','Ecaba2','Eeuro2','Mmula8','Mmuri3','Mmusc10','Mputo1','Ocuni2','Rnorv6','Sscro11']

    chromConverter = {}
    speciesCoordTotRNA = {}
    for species in toChange:
        for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/'+species+'/'+species+'_assembly_report.txt'):
            splitLine = (line.strip()).split('\t')
            if len(splitLine) > 7:
                if not species in ['Dnove3','Eeuro2','Mputo1','Cporc3','Mmuri3','Mputo1','Ocuni2','Sscro11']:
                    chromConverter[splitLine[6]] = 'chr'+str(splitLine[2])
                elif species == 'Cporc3':
                    chromConverter[splitLine[6]] = splitLine[9]
                elif species in ['Mmuri3','Sscro11']:
                    chromConverter[splitLine[4]] = 'chr'+str(splitLine[2])
                elif species == 'Mputo1':
                    chromConverter[splitLine[4]] = splitLine[4].split('.')[0]
                elif species == 'Ocuni2':
                    if splitLine[0].startswith('chr'):
                        chromConverter[splitLine[6]] = splitLine[0]
                    else:
                        chromConverter[splitLine[6]] = 'chr'+str(splitLine[0])
                else:
                    chromConverter[splitLine[6]] = splitLine[4].split('.')[0]

        for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/FINAL_TEST_RUN/NAME_FIX/'+speciesToBed[species]):
            splitLine = (line.strip()).split('\t')
            for i in range(int(splitLine[1]),int(splitLine[2])+1):
                speciesCoordTotRNA[species+'_'+splitLine[0]+'_'+str(i)] = splitLine[3]

    myOutString = ''
    for line in open('/public/groups/corbettlab/tRNA/classifier/JOEL/FINAL_TEST_RUN/all_species_all_trnas.txt'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'species' and partLower(fullToAbb[splitLine[0]]) in toChange:
            myAbb = partLower(fullToAbb[splitLine[0]])
            if myAbb+'_'+chromConverter[splitLine[2]]+'_'+str(int(splitLine[3])+15) in speciesCoordTotRNA:
                mytRNA = speciesCoordTotRNA[myAbb+'_'+chromConverter[splitLine[2]]+'_'+str(int(splitLine[3])+15)]
                splitLine[1] = mytRNA
                splitLine[2] = chromConverter[splitLine[2]]
            else:
                if splitLine[1] == 'KV888377.1.tRNA1-LysTTT':
                    splitLine[1] == 'tRNA-Lys-TTT-1-7'
                    splitLine[2] = chromConverter[splitLine[2]]
                else:
                    print(myAbb+'_'+chromConverter[splitLine[2]]+'_'+str(int(splitLine[3])+15))
        myOutString += joiner(splitLine)+'\n'
    open('reformat_all_trnas.txt','w').write(myOutString)




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





















