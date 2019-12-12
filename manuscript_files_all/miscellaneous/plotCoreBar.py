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

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import matplotlib.image as mpimg

##########################
##### MAIN FUNCTIONS #####
##########################

def plotHumanBar():
    lenToFreq = [0.0]*30
    for line in open('out37.txt'):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 0 and not splitLine[0] == 'Hsapi38':
            lenToFreq[len(splitLine)] += 1

    cumulLenToFreq = [0.0]*30
    runningTotal = 0
    for i in range(2,30):
        cumulLenToFreq[i] = runningTotal + lenToFreq[i]
        runningTotal += lenToFreq[i]


    AAToCoreCount = {}
    for line in open('out37.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0].startswith('Hsapi38-'):
            myAA = splitLine[0].split('-')[2]
            if not myAA in AAToCoreCount:
                AAToCoreCount[myAA] = 0
            AAToCoreCount[myAA] += 1
    AAToDenomCount = {}
    for line in open('FINAL_ORTHOLOG_SET.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[13].startswith('Hsapi38-'):
            myAA = splitLine[13].split('-')[2]
            if not myAA in AAToDenomCount:
                AAToDenomCount[myAA] = 0
            AAToDenomCount[myAA] += 1

    myAAs = []
    myCore = []
    myDenoms = []
    for AA in sorted(AAToDenomCount.keys()):
        if not AA in AAToCoreCount:
            AAToCoreCount[AA] = 0
        myAAs.append(AA)
        myCore.append(AAToCoreCount[AA])
        myDenoms.append(AAToDenomCount[AA])
    print(myAAs)
    print(myCore)
    print(myDenoms)

    fig_width = 13
    fig_height = 10
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.85
    panel_height = 0.75/2

    panel_total_height = (panel_height)*2
    panel_total_width = (panel_width)
    extra_y_space = 1 - panel_total_height
    extra_x_space = 1 - panel_total_width
    above_below = extra_y_space/3
    left_right = extra_x_space/2

    panel1 = plt.axes([left_right, above_below*2+panel_height, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([left_right, above_below, panel_width, panel_height], frameon=True)
    

    panel1.bar(range(0,len(myDenoms)),myDenoms,alpha=1.0,color='gray',label='Human High-Confidence tRNA Genes')
    panel1.bar(range(0,len(myDenoms)),myCore,alpha=1.0,color='black',label='Core tRNA Genes')
    panel1.set_xlim([-1,len(myAAs)])
    panel1.set_xlabel("Amino Acid", fontsize=15)
    panel1.set_ylabel("tRNA Genes", fontsize=15)

    panel2.bar(range(2,30),lenToFreq[2:],edgecolor='black', linewidth=1)
    panel2.set_xlabel("Number of Species with a Syntenic Ortholog", fontsize=15)
    panel2.set_ylabel("Number of Core tRNA Gene Ortholog Sets", fontsize=15)

    panel1.set_xticks(range(0,len(myAAs)))
    panel1.set_xticklabels(myAAs)
    panel1.legend(loc=1,fontsize=15)

    panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=15)
    panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=15)
    plt.savefig('all_human_new.pdf', dpi=700)
    plt.close()


def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return ' '.join(newList)

def main():
    plotHumanBar()

if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit