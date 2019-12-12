#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# separateeCMFiles.py

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

def atacHist():
    cowtRNAToPred = {}
    cowtRNAToChIP = {}
    for f in ['predictions_counts']:
        for line in open('Btaur8_'+f+'.txt'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                cowtRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in cowtRNAToChIP or float(splitLine[2]) > cowtRNAToChIP[splitLine[0]]:
                    cowtRNAToChIP[splitLine[0]] = float(splitLine[2])
    
    goattRNAToPred = {}
    goattRNAToChIP = {}
    for f in ['predictions_counts']:
        for line in open('Chirc1_'+f+'.txt'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                goattRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in goattRNAToChIP or float(splitLine[2]) > goattRNAToChIP[splitLine[0]]:
                    goattRNAToChIP[splitLine[0]] = float(splitLine[2])

    pigtRNAToPred = {}
    pigtRNAToChIP = {}
    for f in ['predictions_counts']:
        for line in open('Sscro11_'+f+'.txt'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                pigtRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in pigtRNAToChIP or float(splitLine[2]) > pigtRNAToChIP[splitLine[0]]:
                    pigtRNAToChIP[splitLine[0]] = float(splitLine[2])

    cowHist1 = []
    cowHist2 = []
    for k in cowtRNAToPred:
        if cowtRNAToChIP[k] == 0:
            if cowtRNAToPred[k] < 0:
                cowHist1.append(getNegIndex(cowtRNAToPred[k]))
            else:
                cowHist1.append(getPosIndex(cowtRNAToPred[k]))
        else:
            if cowtRNAToPred[k] < 0:
                cowHist2.append(getNegIndex(cowtRNAToPred[k]))
            else:
                cowHist2.append(getPosIndex(cowtRNAToPred[k]))

    goatHist1 = []
    goatHist2 = []
    for k in goattRNAToPred:
        if goattRNAToChIP[k] == 0:
            if goattRNAToPred[k] < 0:
                goatHist1.append(getNegIndex(goattRNAToPred[k]))
            else:
                goatHist1.append(getPosIndex(goattRNAToPred[k]))
        else:
            if goattRNAToPred[k] < 0:
                goatHist2.append(getNegIndex(goattRNAToPred[k]))
            else:
                goatHist2.append(getPosIndex(goattRNAToPred[k]))

    pigHist1 = []
    pigHist2 = []
    for k in pigtRNAToPred:
        if pigtRNAToChIP[k] == 0:
            if pigtRNAToPred[k] < 0:
                pigHist1.append(getNegIndex(pigtRNAToPred[k]))
            else:
                pigHist1.append(getPosIndex(pigtRNAToPred[k]))
        else:
            if pigtRNAToPred[k] < 0:
                pigHist2.append(getNegIndex(pigtRNAToPred[k]))
            else:
                pigHist2.append(getPosIndex(pigtRNAToPred[k]))

    fig_width = 16
    fig_height = 14
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.75/2
    panel_height = 0.7/3

    panel_total_height = (panel_height*3)
    panel_total_width = (panel_width*2)
    extra_y_space = 1 - panel_total_height
    extra_x_space = 1 - panel_total_width
    above_below = extra_y_space/4
    left_right = extra_x_space/3

    print(len(cowHist2))
    print(len(goatHist2))
    print(len(pigHist2))

    print(len(cowHist1)+len(cowHist2))
    print(len(goatHist1)+len(goatHist2))
    print(len(pigHist1)+len(pigHist2))


    panel1 = plt.axes([left_right, panel_height*2+above_below*3, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([panel_width+left_right*2, panel_height*2+above_below*3, panel_width, panel_height], frameon=True)
    panel3 = plt.axes([left_right, panel_height+above_below*2, panel_width, panel_height], frameon=True)
    panel4 = plt.axes([panel_width+left_right*2, panel_height+above_below*2, panel_width, panel_height], frameon=True)
    panel5 = plt.axes([left_right, above_below, panel_width, panel_height], frameon=True)
    panel6 = plt.axes([panel_width+left_right*2, above_below, panel_width, panel_height], frameon=True)

    panel1.hist(cowHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel2.hist(cowHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel3.hist(goatHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel4.hist(goatHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel5.hist(pigHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel6.hist(pigHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)

    panel1.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel1.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel2.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel2.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel3.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel3.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel4.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel4.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel5.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel5.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel6.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel6.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])

    panel1.plot([0,0],[-1,500],color='black',linewidth=1)
    panel2.plot([0,0],[-1,500],color='black',linewidth=1)
    panel3.plot([0,0],[-1,500],color='black',linewidth=1)
    panel4.plot([0,0],[-1,500],color='black',linewidth=1)
    panel5.plot([0,0],[-1,500],color='black',linewidth=1)
    panel6.plot([0,0],[-1,500],color='black',linewidth=1)

    panel1.set_xlim([-5.7,5.7])
    panel2.set_xlim([-5.7,5.7])
    panel3.set_xlim([-5.7,5.7])
    panel4.set_xlim([-5.7,5.7])
    panel5.set_xlim([-5.7,5.7])
    panel6.set_xlim([-5.7,5.7])

    panel1.set_ylim([0,225])
    panel2.set_ylim([0,225])
    panel3.set_ylim([0,225])
    panel4.set_ylim([0,225])
    panel5.set_ylim([0,225])
    panel6.set_ylim([0,225])

    panel1.set_ylabel("tRNAs With 0 ATAC-seq peaks", fontsize=16)
    panel2.set_ylabel("tRNAs With > 0 ATAC-seq peaks", fontsize=16)
    panel3.set_ylabel("tRNAs With 0 ATAC-seq peaks", fontsize=16)
    panel4.set_ylabel("tRNAs With > 0 ATAC-seq peaks", fontsize=16)
    panel5.set_ylabel("tRNAs With 0 ATAC-seq peaks", fontsize=16)
    panel6.set_ylabel("tRNAs With > 0 ATAC-seq peaks", fontsize=16)

    panel1.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel1.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel2.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel2.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel3.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel3.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel4.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel4.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel5.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel5.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel6.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel6.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)

    panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel3.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel4.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel5.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel6.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)

    panel1.text(-6.4, 233, 'A', ha='center', va='bottom', fontsize=25)
    panel2.text(-6.4, 233, 'B', ha='center', va='bottom', fontsize=25)
    panel3.text(-6.4, 233, 'C', ha='center', va='bottom', fontsize=25)
    panel4.text(-6.4, 233, 'D', ha='center', va='bottom', fontsize=25)
    panel5.text(-6.4, 233, 'E', ha='center', va='bottom', fontsize=25)
    panel6.text(-6.4, 233, 'F', ha='center', va='bottom', fontsize=25)
    panel1.text(-3.5, 190, 'Cow', ha='center', va='bottom', fontsize=20)
    panel2.text(-3.5, 190, 'Cow', ha='center', va='bottom', fontsize=20)
    panel3.text(-3.5, 190, 'Goat', ha='center', va='bottom', fontsize=20)
    panel4.text(-3.5, 190, 'Goat', ha='center', va='bottom', fontsize=20)
    panel5.text(-3.5, 190, 'Pig', ha='center', va='bottom', fontsize=20)
    panel6.text(-3.5, 190, 'Pig', ha='center', va='bottom', fontsize=20)
    plt.savefig('S7.pdf', dpi=1100)
    plt.close()




def getNegIndex(myPred):
    if myPred > -0.9:
        if myPred > -0.8:
            if myPred > -0.7:
                if myPred > -0.6:
                    return(-1)
                else:
                    return(-2)
            else:
                return(-3)
        else:
            return(-4)
    else:
        return(-5)

def getPosIndex(myPred):
    if myPred < 0.9:
        if myPred < 0.8:
            if myPred < 0.7:
                if myPred < 0.6:
                    return(1)
                else:
                    return(2)
            else:
                return(3)
        else:
            return(4)
    else:
        return(5)





def main():
    reviewerSup()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit











