#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# separateCMFiles.py

import sys
import scipy
import matplotlib
import os
import time

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg

##########################
##### MAIN FUNCTIONS #####
##########################

def plotMouseNew():
    tRNAToNumActive = {}
    for line in open('mmchromatin.txt'):
        splitLine = (line.strip()).split()
        if not splitLine[0] == 'tRNAname':
            tRNAToNumActive[splitLine[0]] = 0
            for i in range(1,len(splitLine)):
                if ',' in splitLine[i]:
                    myStates = splitLine[i].split(',')
                    if '5' in myStates or '7' in myStates or 5 in myStates or 7 in myStates:
                        tRNAToNumActive[splitLine[0]] += 1
                else:
                    if str(splitLine[i]) in ['5','7']:
                        tRNAToNumActive[splitLine[0]] += 1

    tRNAToChrName = {}
    chrNameTotRNA = {}
    for line in open('mm10_activityCoords.txt'):
        splitLine = (line.strip()).split()
        if not splitLine[0] == 'halName':
            tRNAToChrName[splitLine[1]] = splitLine[0]
            chrNameTotRNA[splitLine[0]] = splitLine[1]

    tRNAToNCName = {}
    for line in open('Mmusc10NameMap.txt'):
        splitLine = (line.strip()).split()
        if splitLine[1] in chrNameTotRNA:
            mytRNAName = chrNameTotRNA[splitLine[1]]
            tRNAToNCName[mytRNAName] = splitLine[0]
            #print(mytRNAName, splitLine[0])

    NCToActivity = {}
    for tRNA in tRNAToNCName:
        myActivity = tRNAToNumActive[tRNA]
        myNC = tRNAToNCName[tRNA]
        NCToActivity[myNC] = myActivity

    NCToCpGD = {}
    NCToCpGU = {}
    for line in open('mousetRNADataCpGFixed.tsv'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            NCToCpGD[splitLine[0]] = float(splitLine[3])
            NCToCpGU[splitLine[0]] = float(splitLine[4])

    NCToColor = {}
    for line in open('Mmusc10tRNAClassificationsNewFixedNoSegDups.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[2] == 'active':
            NCToColor[splitLine[0]] = 'red'
        else:
            NCToColor[splitLine[0]] = 'blue'

    myActivities = []
    myCpGDs = []
    myCpGUs = []
    myColors = []
    for tRNA in sorted(NCToCpGU.keys()):
        if tRNA in NCToActivity:
            myActivities.append(NCToActivity[tRNA])
            myCpGDs.append(NCToCpGD[tRNA])
            myCpGUs.append(NCToCpGU[tRNA])
            myColors.append(NCToColor[tRNA])

    myActivitiesNoZero = []
    myCpGDsNoZero = []
    myCpGUsNoZero = []
    for tRNA in sorted(NCToCpGU.keys()):
        #if tRNA in NCToActivity and NCToActivity[tRNA] > 0:
        myActivitiesNoZero.append(NCToActivity[tRNA])
        myCpGDsNoZero.append(NCToCpGD[tRNA])
        myCpGUsNoZero.append(NCToCpGU[tRNA])

    print(scipy.stats.spearmanr(myActivitiesNoZero, myCpGDsNoZero))
    print(scipy.stats.spearmanr(myActivitiesNoZero, myCpGUsNoZero))




    fig_width = 9#8
    fig_height = 9
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.34#0.75
    panel_height = 0.62

    panel_total_height = (panel_height*1)
    panel_total_width = (panel_width*2)
    extra_y_space = 1 - panel_total_height
    extra_x_space = 1 - panel_total_width
    above_below = extra_y_space/2
    left_right = extra_x_space/3

    panel1 = plt.axes([left_right, above_below, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([panel_width+left_right*2, above_below, panel_width, panel_height], frameon=True)

    panel1.scatter(myActivities, myCpGDs, facecolors='none', alpha=0.3, edgecolors=myColors, s=27)
    panel2.scatter(myActivities, myCpGUs, facecolors='none', alpha=0.3, edgecolors=myColors, s=27)
    panel1.set_xlabel('Actively Transcribed Tissues', fontsize=12)
    panel2.set_xlabel('Actively Transcribed Tissues', fontsize=12)
    panel1.set_ylabel('CpG Density Across tRNA Locus', fontsize=12)
    panel2.set_ylabel('CpG Islands Score Upstream of tRNA Gene', fontsize=12)
    panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=12)
    panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=12)
    panel1.text(-1, 0.125, 'A', ha='center', va='bottom', fontsize=26)
    panel2.text(-1,1.25, 'B', ha='center', va='bottom', fontsize=26)
    plt.savefig('reviewerMouseGradient.pdf', dpi=1100)
    plt.close()




def main():
    plotMouseNew()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit













