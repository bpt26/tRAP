#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 10/16/2018
# oneToOneMap.py

import sys
import os
import time
import random
import numpy
import gzip
import math

# build a one-to-one ortholog map across primates in cactus graph from intermediate data
# take the left-most aligning tRNA
# each tRNA must have a reciprocal bond to ONLY ONE OTHER TRNA PER SPECIES

##########################
##### MAIN FUNCTIONS #####
##########################

def oneToOneMap():

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
    for k in sorted(abbrevToFull.keys()):
        if '_' in k:
            mySpeciesNames.append(k)
        else:
            mySpeciesAbbrevs.append(k)

    speciesToDicts = {}
    speciesTotRNAs = {}
    for species1 in mySpeciesAbbrevs:
        speciesTotRNAs[species1] = []
        for species2 in mySpeciesAbbrevs:
            speciesToDicts[species1+'-'+species2] = {}

    mytRNAs = {}

    for line in open('cactustRNAOrthologMap.txt'):
        splitLine = (line.strip()).split('\t')
        tRNA1 = splitLine[0]
        species1 = tRNA1.split('-')[0]
        (speciesTotRNAs[species1]).append(tRNA1)
        if len(splitLine) > 2:
            i = 1
            thisLinetRNAs = [tRNA1]
            while i < len(splitLine):
                tRNA2 = splitLine[i]
                species2 = tRNA2.split('-')[0]
                if float(splitLine[i+1]) > 50:
                    thisLinetRNAs.append(tRNA2)
                i += 2

            for tRNA1 in thisLinetRNAs:
                species1 = tRNA1.split('-')[0]
                for tRNA2 in thisLinetRNAs:
                    species2 = tRNA2.split('-')[0]
                    if not tRNA1 == tRNA2:
                        if not tRNA1 in (speciesToDicts[species1+'-'+species2]) and not tRNA2 in (speciesToDicts[species2+'-'+species1]):
                            (speciesToDicts[species1+'-'+species2])[tRNA1] = tRNA2
                            (speciesToDicts[species2+'-'+species1])[tRNA2] = tRNA1

    printedtRNAs = []
    myOut = 'oneToOneMapAllSpecies.txt'
    open(myOut, 'w').write(joiner(mySpeciesAbbrevs)+'\n')
    for species1 in mySpeciesAbbrevs:
        for tRNA in speciesTotRNAs[species1]:
            if not tRNA in printedtRNAs:
                toPrint = []
                for species2 in mySpeciesAbbrevs:
                    if species1 == species2:
                        toPrint.append(tRNA)
                    elif tRNA in speciesToDicts[species1+'-'+species2] and not (speciesToDicts[species1+'-'+species2])[tRNA] in printedtRNAs:
                        toPrint.append((speciesToDicts[species1+'-'+species2])[tRNA])
                    else:
                        toPrint.append('N/A')
                for i in range(0,len(toPrint)):
                    if toPrint[i] != 'N/A' and toPrint[i] in printedtRNAs:
                        toPrint[i] == 'N/A'
                    printedtRNAs.append(toPrint[i])
                open(myOut, 'a').write(joiner(toPrint)+'\n')


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
    oneToOneMap()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit





