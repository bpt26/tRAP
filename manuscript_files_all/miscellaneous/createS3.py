#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 7/31/2017
# messerBootstrap.py

import sys
import os
import time
import random
import numpy
import gzip
import math

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg


#########################
##### MAIN FUNCTION #####
#########################

def plotFeatureDistributions():
    myActives = []
    myInactives = []
    indexToCat = {}
    catToIndex = {}
    activeCatToVectors = {}
    inactiveCatToVectors = {}
    for line in open('humanCpGTrainingSetFixed.tsv'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'tRNA':
            myHeader = splitLine
            for i in range(1,len(splitLine)-1):
                indexToCat[i] = splitLine[i]
                catToIndex[splitLine[i]] = i
                activeCatToVectors[splitLine[i]] = []
                inactiveCatToVectors[splitLine[i]] = []
        else:
            if str(splitLine[-1]) == '0':
                for i in range(1,len(splitLine)-1):
                    if not splitLine[i] == '?':
                        inactiveCatToVectors[indexToCat[i]].append(float(splitLine[i]))
            else:
                for i in range(1,len(splitLine)-1):
                    if not splitLine[i] == '?':
                        activeCatToVectors[indexToCat[i]].append(float(splitLine[i]))

    labelsDict = {}
    labelsDict['tRNAPhyloPAvg'] = 'Average phyloP Score in tRNA Gene Sequence'
    labelsDict['5PhyloPAvg'] = "Average phyloP Score in 5' Flanking Region"
    labelsDict['CpGDensity'] = 'CpG Density Across tRNA Locus'
    labelsDict['ObsExpUp'] = 'CpG Islands Score Upstream of tRNA Gene'
    labelsDict['GenBit'] = 'tRNAscan-SE General Bit Score'
    labelsDict['tRNA10kb'] = 'tRNA Genes within 10 Kilobases'
    labelsDict['Prot75kb'] = 'Exons within 75 Kilobases'
    labelsDict['TTTT'] = 'Distance to TTTT Transcription Termination Sequence'
    labelsDict['Codon'] = 'Total Number of tRNA Genes With Identical Anticodon'
    labelsDict['MFE'] = 'Minimum Free Energy of Canonical tRNA Structure'



    fig_width = 9*2
    fig_height = 2.8*5
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.34/2
    panel_height = 0.57/5
    panelCounter = 0

    for k in ['CpGDensity','ObsExpUp','5PhyloPAvg','tRNA10kb','Prot75kb','Codon','MFE','GenBit','tRNAPhyloPAvg','TTTT'][::-1]:
        panelCounter += 1
        panel_total_height = (panel_height*5)
        panel_total_width = (panel_width*4)
        extra_y_space = 1 - panel_total_height
        extra_x_space = 1 - panel_total_width
        above_below = extra_y_space/6
        left_right = extra_x_space/5
        myMin = min(min(activeCatToVectors[k]),min(inactiveCatToVectors[k]))
        myMax = max(max(activeCatToVectors[k]),max(inactiveCatToVectors[k]))
        xlim0 = myMin - (0.1*abs(myMax-myMin))
        if xlim0 == 0.0 and '10kb' in k:
            xlim0 = -1
        xlim1 = myMax + (0.1*abs(myMax-myMin))
        xlim = [xlim0,xlim1]
        if panelCounter < 6:
            panel1 = plt.axes([left_right, ((panelCounter*(above_below))+((panelCounter-1)*(panel_height))), panel_width, panel_height], frameon=True)
            panel2 = plt.axes([panel_width+left_right*2, ((panelCounter*(above_below))+((panelCounter-1)*(panel_height))), panel_width, panel_height], frameon=True)
        else:
            tempCounter = panelCounter - 5
            panel1 = plt.axes([panel_width*2+left_right*3, ((tempCounter*(above_below))+((tempCounter-1)*(panel_height))), panel_width, panel_height], frameon=True)
            panel2 = plt.axes([panel_width*3+left_right*4, ((tempCounter*(above_below))+((tempCounter-1)*(panel_height))), panel_width, panel_height], frameon=True)



        if 'tRNA10kb' in k:
            panel1.hist(activeCatToVectors[k], bins=np.arange(max(activeCatToVectors[k])+1)-0.5, edgecolor='black', linewidth=1)
            panel2.hist(inactiveCatToVectors[k], bins=np.arange(max(inactiveCatToVectors[k])+1)-0.5, edgecolor='black', linewidth=1)
        elif 'TTTT' in k:
            print(sorted(activeCatToVectors[k]))
            panel1.hist(activeCatToVectors[k], bins=np.arange(0,360,15)-0.5, edgecolor='black', linewidth=1)
            panel2.hist(inactiveCatToVectors[k], bins=np.arange(0,360,15)-0.5, edgecolor='black', linewidth=1)
        else:
            panel1.hist(activeCatToVectors[k], edgecolor='black', linewidth=1)
            panel2.hist(inactiveCatToVectors[k], edgecolor='black', linewidth=1)
        panel1.set_ylim([0,130])
        panel2.set_ylim([0,130])
        
        panel1.set_xlim(xlim)
        panel2.set_xlim(xlim)
        #panel1.set_xlabel(labelsDict[k], fontsize=12)
        panel1.set_ylabel("Active tRNA Genes", fontsize=12)
        #panel2.set_xlabel(labelsDict[k], fontsize=12)
        panel2.set_ylabel("Inactive tRNA Genes", fontsize=12)
        if k == 'TTTT':
            panel1.set_yticks([0,50,100,130])
            panel1.set_yticklabels(['0','50','100','240'])

        panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=12)
        panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=12)
        if labelsDict[k] not in ["Average phyloP Score in 5' Flanking Region", "Minimum Free Energy of Canonical tRNA Secondary Structure"]:
            panel1.text((xlim[1]+(0.13*xlim[1])), 146, labelsDict[k], ha='center', va='bottom', fontsize=18)
        elif labelsDict[k] == "Average phyloP Score in 5' Flanking Region":
            panel1.text((xlim[1]+(1)), 146, labelsDict[k], ha='center', va='bottom', fontsize=18)
        else:
            panel1.text((xlim[1]+(2.4)), 146, labelsDict[k], ha='center', va='bottom', fontsize=18)
    plt.savefig('active_inactive_ALL_NEW.pdf', dpi=300)
    plt.close()


def joiner(entry):
    """
    Helper function to print lists in 
    tab separated format.
    """
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def getMaxHist(myList):
    myDict = {}
    for k in myList:
        if not k in myDict:
            myDict[k] = 0
        myDict[k] += 1
    myMax = 0
    for k in myDict:
        if myDict[k] > myMax:
            myMax = myDict[k]
    return(myMax)



def main():
    plotFeatureDistributions()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit

