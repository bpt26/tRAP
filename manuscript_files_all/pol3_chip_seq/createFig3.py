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

def getAccuracy(myReads, myPred, myThreshold):
    correct = 0
    incorrect = 0
    for i in range(0,len(myReads)):
        if myReads[i] > myThreshold:
            if myPred[i] > 0.6:
                correct += 1
            else:
                incorrect += 1
        else:
            if myPred[i] > 0.6:
                incorrect += 1
            else:
                correct += 1
    return(float(correct)/float(correct+incorrect)*100)



def getMeanAbove100(myList):
    myReturn = []
    for k in myList:
        if k > 85:
            myReturn.append(k)
    return(np.mean(myReturn))

    
def myLog(myList):
    myReturn = []
    for k in myList:
        if k == 0:
            myReturn.append(0)
        else:
            myReturn.append(np.log10(k))
    return(myReturn)

def floatify(myList):
    myReturn = []
    for k in myList:
        myReturn.append(float(k))
    return(myReturn)


def getAverage(myList, myIndices):
    myReturn = 0.0
    for i in myIndices:
        myReturn += float(myList[i])
    return(myReturn/float(len(myIndices)))

def getMedian(myList, myIndices):
    myReturn = []
    for i in myIndices:
        myReturn.append(float(myList[i]))
    return(sorted(myReturn)[1])

def comparetRAX():
    doubleCopies = {}
    for line in open('mm10-tRNAs.bed'):
        splitLine = (line.strip()).split('\t')
        if str(splitLine[3].split('-')[-1]) == '2':
            doubleCopies['-'.join(splitLine[3].split('-')[:-1])] = True

    # print(len(doubleCopies))
    # print(doubleCopies.keys())

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

    NCToVector = {}
    counter = []
    for line in open('mouseDM-tRNA-seqAll-normalizedreadcounts.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] == 'tRNA':
            myHeader = splitLine[1:]
        else:
            if 'wholecounts' in splitLine[0]:
                mytRNA = splitLine[0].split('_')[0]
                if not mytRNA in doubleCopies:
                    if mytRNA+'-1' in tRNAToNCName:
                        counter.append(mytRNA)
                        if mytRNA == 'tRNA-Arg-TCT-4':
                            print(tRNAToNCName[mytRNA+'-1'], 'Y')
                        NCToVector[tRNAToNCName[mytRNA+'-1']] = floatify(splitLine[1:])

    #print('\n'.join(counter))

    tRNAToPred = {}
    for line in open('Mmusc10tRNAClassificationsNew.txt'):
        splitLine = (line.strip()).split('\t')
        if splitLine[0] in NCToVector:
            tRNAToPred[splitLine[0]] = float(splitLine[1])

    tRNAToLiverChip = {}
    for line in open('mouse_liver_1.bed'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            tRNAToLiverChip[splitLine[0]] = float(splitLine[3])
    for line in open('mouse_liver_2.bed'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            tRNAToLiverChip[splitLine[0]] = (tRNAToLiverChip[splitLine[0]]+float(splitLine[3]))/2.0

    for k in tRNAToNCName.keys():
        print(k+'\t'+tRNAToNCName[k])


    myPred = []
    brainPlus = []
    brainMinus = []
    myColors = []
    myLiverChip = []
    liverPlus = []
    liverMinus = []
    liverBrainPlus = []
    print(myHeader)
    for t in tRNAToPred:
        myPred.append(tRNAToPred[t])
        myLiverChip.append(tRNAToLiverChip[t])
        brainPlus.append((getAverage(NCToVector[t], [0,1,2])))
        #brainMinus.append(getAverage(NCToVector[t], [19,20,21]))
        liverPlus.append(getMedian(NCToVector[t], [6,7,8]))
        #liverMinus.append(getAverage(NCToVector[t], [25,26,27]))
        liverBrainPlus.append(max(getAverage(NCToVector[t], [0,1,2]), getAverage(NCToVector[t], [6,7,8])))
        print(t)
        myColors.append('black')
            
    print(len(myPred))
    liverBrainPlusLog = myLog(liverBrainPlus)

    print(str(getAccuracy(liverBrainPlus, myPred, 80.0))+'%')
    print(str(getAccuracy(liverBrainPlus, myPred, 20.0))+'%')
    print(stats.spearmanr(myLiverChip,liverPlus))

    return(liverPlus, brainPlus, myPred, myColors)


def createPlots(tRNANamesBed, tRNABed, tRNAPredictions, myNameMap, assemblyReport, myLiftOver='', activityLevels =''):
    mytRNAs = []
    startTotRNA = {}
    tRNAToCounts = {}
    tRNAToConf = {}
    tRNAToPredict = {}

    chromConverter = {}
    for line in open(assemblyReport):
        splitLine = (line.strip()).split('\t')
        if len(splitLine) > 8:
            chromConverter[splitLine[6]] = splitLine[9]

    for line in open(tRNANamesBed):
        splitLine = (line.strip()).split('\t')
        mytRNA = splitLine[3]
        mytRNAs.append(mytRNA)
        myStartCoord = int(splitLine[1])
        if 'human' in tRNABed.lower():
            myChrom = splitLine[0]
        else:
            if (splitLine[3]).split('.tRNA')[0] in chromConverter:
                myChrom = chromConverter[(splitLine[3]).split('.tRNA')[0]]
            else:
                myChrom = (splitLine[3]).split('.tRNA')[0]
        for m in range(myStartCoord-10, myStartCoord+10):
            myStart = myChrom+':'+str(m)
            startTotRNA[myStart] = mytRNA

    oldNameToStart = {}
    if len(myLiftOver) > 0: 
        for line in open(myLiftOver):
            splitLine = (line.strip()).split('\t')
            myOldName = str(splitLine[3])
            oldNameToStart[myOldName] = str(splitLine[0])+':'+str(splitLine[1])

    nameMap = {}
    for line in open(myNameMap):
        splitLine = (line.strip()).split('\t')
        nameMap[splitLine[0]] = splitLine[1]
        nameMap[splitLine[1]] = splitLine[0]

    for line in open(tRNABed):
        splitLine = (line.strip()).split('\t')
        if len(oldNameToStart) > 1:
            if splitLine[0] in oldNameToStart:
                myStart = oldNameToStart[splitLine[0]]
                if myStart in startTotRNA:
                    mytRNA = startTotRNA[myStart]
                    #print(mytRNA)
                    #mytRNA = nameMap[mytRNA]
                    #print(mytRNA)
                    tRNAToCounts[mytRNA] = int(splitLine[5])
                    tRNAToCounts[nameMap[mytRNA]] = int(splitLine[5])
                # else:
                #     print(myStart)
                #     print(splitLine[0])
        else:
            myStart = str(splitLine[1])+':'+str(splitLine[2])
            if myStart in startTotRNA:
                mytRNA = startTotRNA[myStart]
                #mytRNA = nameMap[mytRNA]
                #print(mytRNA)
                tRNAToCounts[mytRNA] = int(splitLine[5])
            # else:
            #     print(myStart)
            #     print(splitLine[0])
    #print(tRNAToCounts)
    #print(startTotRNA)
    
    for line in open(tRNAPredictions):
        splitLine = (line.strip()).split('\t')
        mytRNA = str(splitLine[0])
        myPredict = str(splitLine[2])
        tRNAToConf[mytRNA] = float(splitLine[1])
        #if not 'mouse' in tRNABed.lower():
        if 'inactive' in myPredict:
            tRNAToPredict[mytRNA] = 'black'
        else:
            tRNAToPredict[mytRNA] = 'black'

    if len(activityLevels) > 0:
        for line in open(activityLevels):
            splitLine = (line.strip()).split('\t')
            if splitLine[1] == 'inactive':
                tRNAToPredict[splitLine[0]] = 'black'#'coral'
            else:
                tRNAToPredict[splitLine[0]] = 'black'#'purple'

    myCounts = []
    myConfs = []
    myPredicts = []
    myOut = tRNABed.split('.')[0]+'.bed'
    open(myOut, 'w').write('tRNA\tconf\tpredict\tcount\n')
    percentNonZero = [0.0,0.0]
    for tRNA in mytRNAs:
        if tRNA in tRNAToConf:
            print(tRNA)
            myConfs.append(tRNAToConf[tRNA])
            myPredicts.append(tRNAToPredict[tRNA])
            if tRNA in tRNAToCounts:
                myCounts.append(tRNAToCounts[tRNA])
                open(myOut, 'a').write(tRNA+'\t'+str(tRNAToConf[tRNA])+'\t'+str(tRNAToPredict[tRNA])+'\t'+str(tRNAToCounts[tRNA])+'\n')
                if tRNAToCounts[tRNA] > 0:
                    percentNonZero[1] += 1.0
                    if tRNAToConf[tRNA] > 0:
                        percentNonZero[0] += 1.0
            else:
                myCounts.append(0)
                open(myOut, 'a').write(tRNA+'\t'+str(tRNAToConf[tRNA])+'\t'+str(tRNAToPredict[tRNA])+'\t0\n')
        else:
            tRNA = nameMap[tRNA]
            if tRNA in tRNAToConf:
                #print(tRNA)
                myConfs.append(tRNAToConf[tRNA])
                myPredicts.append(tRNAToPredict[tRNA])
                if tRNA in tRNAToCounts:
                    myCounts.append(tRNAToCounts[tRNA])
                    open(myOut, 'a').write(tRNA+'\t'+str(tRNAToConf[tRNA])+'\t'+str(tRNAToPredict[tRNA])+'\t'+str(tRNAToCounts[tRNA])+'\n')
                    if tRNAToCounts[tRNA] > 0:
                        percentNonZero[1] += 1.0
                        if tRNAToConf[tRNA] > 0:
                            percentNonZero[0] += 1.0
                else:
                    myCounts.append(0)
                    open(myOut, 'a').write(tRNA+'\t'+str(tRNAToConf[tRNA])+'\t'+str(tRNAToPredict[tRNA])+'\t0\n')
    print((percentNonZero[0]/percentNonZero[1])*100)
    #print(myCounts, myConfs, myPredicts)
    return(myCounts, myConfs, myPredicts)

fig_width = 9.5
fig_height = 7.0
plt.figure(figsize=(fig_width, fig_height))
panel_width = 0.315
panel_height = 0.24

panel_total_height = (panel_height*3)
extra_y_space = 1 - panel_total_height
above_below = extra_y_space/4

panel_total_width = (panel_width*2)
extra_x_space = 1 - panel_total_width
left_right = extra_y_space/3

panel1a = plt.axes([left_right, above_below*3+panel_height*2, panel_width/2, panel_height], frameon=True)
panel1b = plt.axes([left_right+(panel_width/2), above_below*3+panel_height*2, panel_width/2, panel_height], frameon=True)
panel2a = plt.axes([1-panel_width-left_right, above_below*3+panel_height*2, panel_width/2, panel_height], frameon=True)
panel2b = plt.axes([1-panel_width-left_right+(panel_width/2), above_below*3+panel_height*2, panel_width/2, panel_height], frameon=True)

panel3a = plt.axes([left_right, above_below*2+panel_height, panel_width/2, panel_height], frameon=True)
panel3b = plt.axes([left_right+(panel_width/2), above_below*2+panel_height, panel_width/2, panel_height], frameon=True)
panel4a = plt.axes([1-panel_width-left_right, above_below*2+panel_height, panel_width/2, panel_height], frameon=True)
panel4b = plt.axes([1-panel_width-left_right+(panel_width/2), above_below*2+panel_height, panel_width/2, panel_height], frameon=True)

panel5a = plt.axes([left_right, above_below, panel_width/2, panel_height], frameon=True)
panel5b = plt.axes([left_right+(panel_width/2), above_below, panel_width/2, panel_height], frameon=True)
panel6a = plt.axes([1-panel_width-left_right, above_below, panel_width/2, panel_height], frameon=True)
panel6b = plt.axes([1-panel_width-left_right+(panel_width/2), above_below, panel_width/2, panel_height], frameon=True)

myCounts, myConfs, myPredicts = createPlots('Mmusc10_halChroms_hiConftRNAs.bed', 'mouse_liver_1.txt', 'Mmusc10tRNAClassificationsNewFixedNoSegDups.txt', 'Mmusc10NameMap.txt', 'Mmusc10_assembly_report.txt', 'mouseLift.txt', 'mouseActivityInput.txt')
panel1a.set_xlim(-1.05,-0.45)
panel1b.set_xlim(0.45,1.05)
panel1aTicks = [-1, -0.5]
panel1a.set_xticks(panel1aTicks)
panel1a.set_xticklabels(['1', '0.5'])
panel1bTicks = [0.5, 1]
panel1b.set_xticklabels(['0.5', '1'])
panel1b.set_xticks(panel1bTicks)
panel1a.spines['right'].set_visible(False)
panel1b.spines['left'].set_visible(False)
panel1a.tick_params(right='off', labelright='off')
panel1b.tick_params(left='off', labelleft='off')
#panel1a.set_ylabel('Mouse Gene Read Count', fontsize=11)
#panel1a.set_xlabel('Predicted Inactive')
#panel1b.set_xlabel('Predicted Active')
panel1a.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel1b.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel1a.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel1b.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel1a.plot([-0.45, -0.45], [min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))], color=(0,0,0), linewidth=3)
panel1a.text(-1.05, 1.12*max(myCounts), "A", ha='center', va='bottom', color=(0,0,0), fontsize=18)

myCounts, myConfs, myPredicts = createPlots('Mmula8_hiConftRNAs.bed', 'macaque_liver_1.txt', 'Mmula8tRNAClassificationsNewFixedNoSegDups.txt', 'Mmula8NameMap.txt', 'Mmula8_assembly_report.txt', 'rhesusLift.txt')
panel2a.set_xlim(-1.05,-0.45)
panel2b.set_xlim(0.45,1.05)
panel2aTicks = [-1, -0.5]
panel2a.set_xticks(panel2aTicks)
panel2a.set_xticklabels(['1', '0.5'])
panel2bTicks = [0.5, 1]
panel2b.set_xticklabels(['0.5', '1'])
panel2b.set_xticks(panel2bTicks)
panel2a.spines['right'].set_visible(False)
panel2b.spines['left'].set_visible(False)
panel2a.tick_params(right='off', labelright='off')
panel2b.tick_params(left='off', labelleft='off')
#panel2a.set_ylabel('Macaque Gene Read Count', fontsize=11)
#panel2a.set_xlabel('Predicted Inactive')
#panel2b.set_xlabel('Predicted Active')
panel2a.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel2b.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel2a.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel2b.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel2a.plot([-0.45, -0.45], [min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))], color=(0,0,0), linewidth=3)
panel2a.text(-1.05, 1.12*max(myCounts), "B", ha='center', va='bottom', color=(0,0,0), fontsize=18)

myCounts, myConfs, myPredicts = createPlots('Rnorv6_hiConftRNAs.bed', 'rat_liver_1.txt', 'Rnorv6tRNAClassificationsNewFixedNoSegDups.txt', 'Rnorv6NameMap.txt', 'Rnorv6_assembly_report.txt', 'ratLift.txt')
panel3a.set_xlim(-1.05,-0.45)
panel3b.set_xlim(0.45,1.05)
panel3aTicks = [-1, -0.5]
panel3a.set_xticks(panel3aTicks)
panel3a.set_xticklabels(['1', '0.5'])
panel3bTicks = [0.5, 1]
panel3b.set_xticklabels(['0.5', '1'])
panel3b.set_xticks(panel3bTicks)
panel3a.spines['right'].set_visible(False)
panel3b.spines['left'].set_visible(False)
panel3a.tick_params(right='off', labelright='off')
panel3b.tick_params(left='off', labelleft='off')
#panel3a.set_ylabel('Rat Gene Read Count', fontsize=11)
#panel3a.set_xlabel('Predicted Inactive')
#panel3b.set_xlabel('Predicted Active')
panel3a.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel3b.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel3a.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel3b.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel3a.plot([-0.45, -0.45], [min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))], color=(0,0,0), linewidth=3)
panel3a.text(-1.05, 1.12*max(myCounts), "C", ha='center', va='bottom', color=(0,0,0), fontsize=18)
#panel3a.text(-0.45, -0.45*max(myCounts), "Confidence", ha='center', va='bottom', color=(0,0,0), fontsize=13)

myCounts, myConfs, myPredicts = createPlots('Cfami3_hiConftRNAs.bed', 'dog_liver_1.txt', 'Cfami3tRNAClassificationsNewFixedNoSegDups.txt', 'Cfami3NameMap.txt', 'Cfami3_assembly_report.txt', 'dogLift.txt')
panel4a.set_xlim(-1.05,-0.45)
panel4b.set_xlim(0.45,1.05)
panel4aTicks = [-1, -0.5]
panel4a.set_xticks(panel4aTicks)
panel4a.set_xticklabels(['1', '0.5'])
panel4bTicks = [0.5, 1]
panel4b.set_xticklabels(['0.5', '1'])
panel4b.set_xticks(panel4bTicks)
panel4a.spines['right'].set_visible(False)
panel4b.spines['left'].set_visible(False)
panel4a.tick_params(right='off', labelright='off')
panel4b.tick_params(left='off', labelleft='off')
#panel4a.set_ylabel('Dog Gene Read Count', fontsize=11)
#panel4a.set_xlabel('Predicted Inactive')
#panel4b.set_xlabel('Predicted Active')
panel4a.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel4b.scatter(myConfs, myCounts, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel4a.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel4b.set_ylim([min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))])
panel4a.plot([-0.45, -0.45], [min(myCounts)-(.1*max(myCounts)), max(myCounts)+(.1*max(myCounts))], color=(0,0,0), linewidth=3)
panel4a.text(-1.05, 1.12*max(myCounts), "D", ha='center', va='bottom', color=(0,0,0), fontsize=18)
#panel4a.text(-0.45, -0.45*max(myCounts), "Confidence", ha='center', va='bottom', color=(0,0,0), fontsize=13)

liverPlus, brainPlus, myConfs, myPredicts = comparetRAX()

panel5a.set_xlim(-1.05,-0.45)
panel5b.set_xlim(0.45,1.05)
panel5aTicks = [-1, -0.5]
panel5a.set_xticks(panel5aTicks)
panel5a.set_xticklabels(['1', '0.5'])
panel5bTicks = [0.5, 1]
panel5b.set_xticklabels(['0.5', '1'])
panel5b.set_xticks(panel5bTicks)
panel5a.spines['right'].set_visible(False)
panel5b.spines['left'].set_visible(False)
panel5a.tick_params(right='off', labelright='off')
panel5b.tick_params(left='off', labelleft='off')
#panel5a.set_ylabel('Liver DM-tRNA-seq Read Count', fontsize=11)
#panel5a.set_xlabel('Predicted Inactive')
#panel5b.set_xlabel('Predicted Active')
panel5a.scatter(myConfs, liverPlus, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel5b.scatter(myConfs, liverPlus, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel5a.set_ylim([0.0-(.1*5000), 5000+(.1*5000)])
panel5b.set_ylim([0.0-(.1*5000), 5000+(.1*5000)])
panel5a.plot([-0.45, -0.45], [0.0-(.1*5000), 5000+(.1*5000)], color=(0,0,0), linewidth=3)
panel5a.text(-1.05, 1.12*5000, "E", ha='center', va='bottom', color=(0,0,0), fontsize=18)
#panel5a.text(-0.45, -0.45*max(myCounts), "Confidence", ha='center', va='bottom', color=(0,0,0), fontsize=13)

panel6a.set_xlim(-1.05,-0.45)
panel6b.set_xlim(0.45,1.05)
panel6aTicks = [-1, -0.5]
panel6a.set_xticks(panel6aTicks)
panel6a.set_xticklabels(['1', '0.5'])
panel6bTicks = [0.5, 1]
panel6b.set_xticklabels(['0.5', '1'])
panel6b.set_xticks(panel6bTicks)
panel6a.spines['right'].set_visible(False)
panel6b.spines['left'].set_visible(False)
panel6a.tick_params(right='off', labelright='off')
panel6b.tick_params(left='off', labelleft='off')
#panel6a.set_ylabel('Liver DM-tRNA-seq Read Count', fontsize=11)
#panel6a.set_xlabel('Predicted Inactive')
#panel6b.set_xlabel('Predicted Active')
panel6a.scatter(myConfs, brainPlus, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
panel6b.scatter(myConfs, brainPlus, alpha=0.3, facecolors='none', edgecolors=myPredicts, s=18)
plt.subplots_adjust(wspace=0)
panel6a.set_ylim([0.0-(.1*5000), 5000+(.1*5000)])
panel6b.set_ylim([0.0-(.1*5000), 5000+(.1*5000)])
panel6a.plot([-0.45, -0.45], [0.0-(.1*5000), 5000+(.1*5000)], color=(0,0,0), linewidth=3)
panel6a.text(-1.05, 1.12*5000, "F", ha='center', va='bottom', color=(0,0,0), fontsize=18)
#panel6a.text(-0.45, -0.45*max(myCounts), "Confidence", ha='center', va='bottom', color=(0,0,0), fontsize=13)


myFileName = 'predictionsFourPanelRFLiverNewWithDM.pdf'
plt.savefig(myFileName, dpi=1100)

