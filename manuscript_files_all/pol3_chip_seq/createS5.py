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

def ChIPHist():
    mousetRNAToPred = {}
    mousetRNAToChIP = {}
    for f in ['liver_1','liver_2','muscle_1','muscle_2','testis_1','testis_2']:
        for line in open('mouse_'+f+'.bed'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                mousetRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in mousetRNAToChIP or int(splitLine[3]) > mousetRNAToChIP[splitLine[0]]:
                    mousetRNAToChIP[splitLine[0]] = int(splitLine[3])
    
    macaquetRNAToPred = {}
    macaquetRNAToChIP = {}
    for f in ['liver_1','liver_2']:
        for line in open('macaque_'+f+'.bed'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                macaquetRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in macaquetRNAToChIP or int(splitLine[3]) > macaquetRNAToChIP[splitLine[0]]:
                    macaquetRNAToChIP[splitLine[0]] = int(splitLine[3])

    rattRNAToPred = {}
    rattRNAToChIP = {}
    for f in ['liver_1','liver_2']:
        for line in open('rat_'+f+'.bed'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                rattRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in rattRNAToChIP or int(splitLine[3]) > rattRNAToChIP[splitLine[0]]:
                    rattRNAToChIP[splitLine[0]] = int(splitLine[3])
    
    dogtRNAToPred = {}
    dogtRNAToChIP = {}
    for f in ['liver_1','liver_2']:
        for line in open('dog_'+f+'.bed'):
            splitLine = (line.strip().split('\t'))
            if not splitLine[0] == 'tRNA':
                dogtRNAToPred[splitLine[0]] = float(splitLine[1])
                if not splitLine[0] in dogtRNAToChIP or int(splitLine[3]) > dogtRNAToChIP[splitLine[0]]:
                    dogtRNAToChIP[splitLine[0]] = int(splitLine[3])

    mouseHist1 = []
    mouseHist2 = []
    mouse1 = 0
    mouse2 = 0
    for k in mousetRNAToPred:
        if mousetRNAToChIP[k] == 0:
            if mousetRNAToPred[k] < 0:
                mouseHist1.append(getNegIndex(mousetRNAToPred[k]))
            else:
                mouseHist1.append(getPosIndex(mousetRNAToPred[k]))
                mouse2 += 1
                
        else:
            mouse1 += 1
            if mousetRNAToPred[k] < 0:
                mouseHist2.append(getNegIndex(mousetRNAToPred[k]))
            else:
                mouseHist2.append(getPosIndex(mousetRNAToPred[k]))

    macaqueHist1 = []
    macaqueHist2 = []
    mac1 = 0
    mac2 = 0
    for k in macaquetRNAToPred:
        if macaquetRNAToChIP[k] == 0:
            if macaquetRNAToPred[k] < 0:
                macaqueHist1.append(getNegIndex(macaquetRNAToPred[k]))
            else:
                macaqueHist1.append(getPosIndex(macaquetRNAToPred[k]))
                mac2 += 1
                
        else:
            mac1 += 1
            if macaquetRNAToPred[k] < 0:
                macaqueHist2.append(getNegIndex(macaquetRNAToPred[k]))
            else:
                macaqueHist2.append(getPosIndex(macaquetRNAToPred[k]))

    ratHist1 = []
    ratHist2 = []
    rat1 = 0
    rat2 = 0
    for k in rattRNAToPred:
        if rattRNAToChIP[k] == 0:
            if rattRNAToPred[k] < 0:
                ratHist1.append(getNegIndex(rattRNAToPred[k]))
            else:
                ratHist1.append(getPosIndex(rattRNAToPred[k]))
                rat2 += 1
                
        else:
            rat1 += 1
            if rattRNAToPred[k] < 0:
                ratHist2.append(getNegIndex(rattRNAToPred[k]))
            else:
                ratHist2.append(getPosIndex(rattRNAToPred[k]))

    dogHist1 = []
    dogHist2 = []
    dog1 = 0
    dog2 = 0
    for k in dogtRNAToPred:
        if dogtRNAToChIP[k] == 0:
            if dogtRNAToPred[k] < 0:
                dogHist1.append(getNegIndex(dogtRNAToPred[k]))
            else:
                dogHist1.append(getPosIndex(dogtRNAToPred[k]))
                dog2 += 1
                
        else:
            dog1 += 1
            if dogtRNAToPred[k] < 0:
                dogHist2.append(getNegIndex(dogtRNAToPred[k]))
            else:
                dogHist2.append(getPosIndex(dogtRNAToPred[k]))

    fig_width = 16
    fig_height = 18
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.75/2
    panel_height = 0.7/4

    panel_total_height = (panel_height*4)
    panel_total_width = (panel_width*2)
    extra_y_space = 1 - panel_total_height
    extra_x_space = 1 - panel_total_width
    above_below = extra_y_space/5
    left_right = extra_x_space/3

    print(mouse1, mouse2, len(mouseHist1), len(mouseHist1)+len(mouseHist2))
    print(mac1, mac2, len(macaqueHist1), len(macaqueHist1)+len(macaqueHist2))
    print(rat1, rat2, len(ratHist1), len(ratHist1)+len(ratHist2))
    print(dog1, dog2, len(dogHist1), len(dogHist1)+len(dogHist2))

    panel1 = plt.axes([left_right, panel_height*3+above_below*4, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([panel_width+left_right*2, panel_height*3+above_below*4, panel_width, panel_height], frameon=True)
    panel3 = plt.axes([left_right, panel_height*2+above_below*3, panel_width, panel_height], frameon=True)
    panel4 = plt.axes([panel_width+left_right*2, panel_height*2+above_below*3, panel_width, panel_height], frameon=True)

    panel5 = plt.axes([left_right, panel_height+above_below*2, panel_width, panel_height], frameon=True)
    panel6 = plt.axes([panel_width+left_right*2, panel_height+above_below*2, panel_width, panel_height], frameon=True)
    panel7 = plt.axes([left_right, above_below, panel_width, panel_height], frameon=True)
    panel8 = plt.axes([panel_width+left_right*2, above_below, panel_width, panel_height], frameon=True)

    panel1.hist(mouseHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel2.hist(mouseHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel3.hist(macaqueHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel4.hist(macaqueHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel5.hist(ratHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel6.hist(ratHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel7.hist(dogHist1,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)
    panel8.hist(dogHist2,edgecolor='black', linewidth=1, bins=np.arange(-5,7)-0.5)

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
    panel7.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel7.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])
    panel8.set_xticks([-5,-4,-3,-2,-1,1,2,3,4,5])
    panel8.set_xticklabels(['>0.9','>0.8','>0.7','>0.6','>0.5','>0.5','>0.6','>0.7','>0.8','>0.9'])

    panel1.plot([0,0],[-1,500],color='black',linewidth=1)
    panel2.plot([0,0],[-1,500],color='black',linewidth=1)
    panel3.plot([0,0],[-1,500],color='black',linewidth=1)
    panel4.plot([0,0],[-1,500],color='black',linewidth=1)
    panel5.plot([0,0],[-1,500],color='black',linewidth=1)
    panel6.plot([0,0],[-1,500],color='black',linewidth=1)
    panel7.plot([0,0],[-1,500],color='black',linewidth=1)
    panel8.plot([0,0],[-1,500],color='black',linewidth=1)

    panel1.set_xlim([-5.7,5.7])
    panel2.set_xlim([-5.7,5.7])
    panel3.set_xlim([-5.7,5.7])
    panel4.set_xlim([-5.7,5.7])
    panel5.set_xlim([-5.7,5.7])
    panel6.set_xlim([-5.7,5.7])
    panel7.set_xlim([-5.7,5.7])
    panel8.set_xlim([-5.7,5.7])

    panel1.set_ylim([0,200])
    panel2.set_ylim([0,200])
    panel3.set_ylim([0,200])
    panel4.set_ylim([0,200])
    panel5.set_ylim([0,200])
    panel6.set_ylim([0,200])
    panel7.set_ylim([0,200])
    panel8.set_ylim([0,200])

    panel1.set_ylabel("tRNAs With 0 ChIP-seq reads", fontsize=16)
    panel2.set_ylabel("tRNAs With > 0 ChIP-seq reads", fontsize=16)
    panel3.set_ylabel("tRNAs With 0 ChIP-seq reads", fontsize=16)
    panel4.set_ylabel("tRNAs With > 0 ChIP-seq reads", fontsize=16)
    panel5.set_ylabel("tRNAs With 0 ChIP-seq reads", fontsize=16)
    panel6.set_ylabel("tRNAs With > 0 ChIP-seq reads", fontsize=16)
    panel7.set_ylabel("tRNAs With 0 ChIP-seq reads", fontsize=16)
    panel8.set_ylabel("tRNAs With > 0 ChIP-seq reads", fontsize=16)

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
    panel7.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel7.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)
    panel8.text(-3, -36, 'Predicted Inactive', ha='center', va='bottom', fontsize=14)
    panel8.text(3, -36, 'Predicted Active', ha='center', va='bottom', fontsize=14)

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
    panel7.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel8.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=14)
    panel1.text(-6.4, 213, 'A', ha='center', va='bottom', fontsize=25)
    panel2.text(-6.4, 213, 'B', ha='center', va='bottom', fontsize=25)
    panel3.text(-6.4, 213, 'C', ha='center', va='bottom', fontsize=25)
    panel4.text(-6.4, 213, 'D', ha='center', va='bottom', fontsize=25)
    panel5.text(-6.4, 213, 'E', ha='center', va='bottom', fontsize=25)
    panel6.text(-6.4, 213, 'F', ha='center', va='bottom', fontsize=25)
    panel7.text(-6.4, 213, 'G', ha='center', va='bottom', fontsize=25)
    panel8.text(-6.4, 213, 'H', ha='center', va='bottom', fontsize=25)
    panel1.text(-3.5, 175, 'Mouse', ha='center', va='bottom', fontsize=20)
    panel2.text(-3.5, 175, 'Mouse', ha='center', va='bottom', fontsize=20)
    panel3.text(-3.5, 175, 'Macaque', ha='center', va='bottom', fontsize=20)
    panel4.text(-3.5, 175, 'Macaque', ha='center', va='bottom', fontsize=20)
    panel5.text(-3.5, 175, 'Rat', ha='center', va='bottom', fontsize=20)
    panel6.text(-3.5, 175, 'Rat', ha='center', va='bottom', fontsize=20)
    panel7.text(-3.5, 175, 'Dog', ha='center', va='bottom', fontsize=20)
    panel8.text(-3.5, 175, 'Dog', ha='center', va='bottom', fontsize=20)
    plt.savefig('S4.pdf', dpi=1100)
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
    ChIPHist()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit











