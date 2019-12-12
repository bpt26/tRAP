#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 11/16/2017
# classifierv4.py

import sys
import scipy
import matplotlib
import sklearn
from sklearn import svm
import numpy as np
import matplotlib.pyplot as plt

from sklearn import decomposition
from sklearn import datasets
from sklearn import tree
from sklearn.ensemble import RandomForestClassifier
from sklearn.datasets import make_classification
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import train_test_split
from sklearn.model_selection import cross_validate
from sklearn.model_selection import cross_val_predict
from sklearn.impute import SimpleImputer
from sklearn.metrics import roc_curve, auc
from sklearn.svm import SVC

def buildClassifier():

    labelsDict = {}
    labelsDict['tRNAPhyloPAvg'] = 'Average PhyloP\nScore in tRNA\nGene Sequence'
    labelsDict['5PhyloPAvg'] = "Average PhyloP\nScore in\n5' Flanking Region"
    labelsDict['CpGOvrPct'] = 'Percentage of CpG\nDinucleotides Across\ntRNA Locus'
    labelsDict['ObsExp'] = 'Observed/Expected\nCpG Islands Score\nAcross tRNA Locus'
    labelsDict['ObsExpUp'] = 'Observed/Expected\nCpG Islands Score\nUpstream of tRNA Gene'
    labelsDict['GenBit'] = 'tRNAscan-SE General\nBit Score'
    labelsDict['tRNA10kb'] = 'tRNA Genes within 10 Kilobases'
    labelsDict['Prot75kb'] = 'Exons within 75 Kilobases'
    labelsDict['TTTT'] = 'Distance to Nearest TTTT\nTranscription Termination\nSequence'
    labelsDict['Codon'] = 'tRNAs Corresponding\nto the Same Codon'
    labelsDict['MFE'] = 'Constrained Minimum\nFree Energy'

    myHumanData = []
    myLabels = []
    myHumanNames = []
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    for line in open('humanCpGTrainingSet.tsv'):
        splitLine = (line.strip()).split('\t')
        if (splitLine[0]) == 'tRNA':
            myHeader = []
            for k in splitLine[1:]:
                myHeader.append(k)
        else:
            myHumanData.append(makeFloat(splitLine[1:-1]))
            myHumanNames.append(splitLine[0])
            if str(splitLine[-1]) in ['active','1']:
                myLabels.append(1)
            elif str(splitLine[-1]) in ['inactive','0']:
                myLabels.append(0)

    imp_mean.fit_transform(myHumanData)
    myHumanDataReplaced = imp_mean.transform(myHumanData)

    myMouseData = []
    myMouseLabels = []
    myMouseNames = []
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    for line in open('mousetRNADataCpG.tsv'):
        splitLine = (line.strip()).split('\t')
        if (splitLine[0]) == 'tRNA':
            myMouseHeader = splitLine
        else:
            myMouseData.append(makeFloat(splitLine[1:-1]))
            myMouseNames.append(splitLine[0])
            if str(splitLine[-1]) in ['active','1']:
                myMouseLabels.append(1)
            elif str(splitLine[-1]) in ['inactive','0']:
                myMouseLabels.append(0)

    imp_mean.fit_transform(myMouseData)
    myMouseDataReplaced = imp_mean.transform(myMouseData)

    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    #print(clf.score(myMouseDataReplaced,myMouseLabels))
    myPredictions = clf.predict(myMouseDataReplaced)
    for i in range(0,len(myPredictions)):
        if not myPredictions[i] == myMouseLabels[i]:
            print(myMouseNames[i], myPredictions[i], myMouseLabels[i])
    cM = confusionMatrix(clf.predict(myMouseDataReplaced),myMouseLabels)
    # for i in range(0,len(myMouseLabels)):
    #     if not myPredictions[i] == myMouseLabels[i]:
    #         print(myMouseNames[i], myPredictions[i], myMouseLabels[i])

    print(cM)
    print(getAccuracy(cM))
    print(getScore(clf.predict_proba(myMouseDataReplaced),myMouseLabels))


    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    myPredictions = clf.predict_proba(myMouseDataReplaced)
    myOutString = ''
    for i in range(0,len(myPredictions)):
        myOutString += myMouseNames[i]+'\t'
        if float(myPredictions[i][0]) > float(myPredictions[i][1]):
            myOutString += '-'+str(myPredictions[i][0])+'\tinactive\n'
        else:
            myOutString += str(myPredictions[i][1])+'\tactive\n'
    open('mousePredictionsNew.txt', 'w').write(myOutString)



    """
    Our final model uses a bag size of 100%, 200 iterations, evaluation of 2 attributes at each node,
    a minimum variance of 1e-4 per split, and a maximum depth of 5 nodes.
    """

    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions1 = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict')
    print(cvPredictions1)
    cM = confusionMatrix(cvPredictions1, myLabels)
    for i in range(0,len(cvPredictions1)):
        if not cvPredictions1[i] == myLabels[i]:
            print(myHumanNames[i], cvPredictions1[i], myLabels[i])
    print(cM)
    print(getAccuracy(cM))
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    #print(cvPredictions)
    print(getScore(cvPredictions, myLabels))


    myOutString = ''
    for i in range(0,len(cvPredictions)):
        myOutString += myHumanNames[i]+'\t'
        if float(cvPredictions[i][0]) > float(cvPredictions[i][1]):
            myOutString += '-'+str(cvPredictions[i][0])+'\tinactive\n'
        else:
            myOutString += str(cvPredictions[i][1])+'\tactive\n'
    open('humanCVPredictions.txt', 'w').write(myOutString)




    
    fig_width = 8
    fig_height = 8
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.63
    panel_height = 0.63
    panel1 = plt.axes([(1.0-panel_width)/2, (1-panel_height)/2, panel_width, panel_height], frameon=True)
    myColors = ['red','orange','yellow','green','cyan','blue','indigo','violet','brown','pink']
    print(clf.feature_importances_)
    # panel1.pie(reorder(clf.feature_importances_,[0,1,6,7,5,2,3,4,8,9]), labels=reorder(myHeader,[0,1,6,7,5,2,3,4,8,9]), colors=myColors)
    # plt.axis('equal')
    # plt.tick_params(bottom='on', labelbottom='on',\
    #            left='on', labelleft='on', \
    #            right='off', labelright='off',\
    #            top='off', labeltop='off', labelsize=12)
    # plt.savefig('featureImportancePieChanged.png', dpi=700)
    # plt.close()


    # cM = confusionMatrix(cvPredictions, np.asarray(myLabels))
    # print(cM)
    # print(getAccuracy(cM))


    fig_width = 21
    fig_height = 7
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.4*(2.0/3.0)
    panel_height = 0.8
    panel_total_height = (panel_height*1)
    extra_y_space = 1 - panel_total_height
    above_below = extra_y_space/2
    panel_total_width = (panel_width*2)
    extra_x_space = 1 - panel_total_width
    left_right = extra_y_space/3
    panel1 = plt.axes([left_right, above_below, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([panel_width+(2*left_right), above_below, panel_width, panel_height], frameon=True)


    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    fpr, tpr, thresholds = roc_curve(myLabels, cvPredictions[:,1])
    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='b', label='Random Forest (AUC = %0.3f)' % (roc_auc))

    clf = LogisticRegression(random_state=49, solver='lbfgs', multi_class='multinomial')
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    fpr, tpr, thresholds = roc_curve(myLabels, cvPredictions[:,1])
    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='r', label='Logistic Regression (AUC = %0.3f)' % (roc_auc))
    print(clf.coef_)

    clf = SVC(probability=True, gamma='auto', kernel='linear')
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')


    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='y', label='Support Vector Machine (AUC = %0.3f)' % (roc_auc))

    panel1.set_xlabel("False Positive Rate", fontsize=18)
    panel1.set_ylabel("True Positive Rate", fontsize=18)

    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    mousePred = clf.predict_proba(myMouseDataReplaced)
    fpr, tpr, thresholds = roc_curve(myMouseLabels, mousePred[:,1])
    roc_auc = auc(fpr, tpr)
    panel2.plot(fpr, tpr, color='b', label='Random Forest (AUC = %0.3f)' % (roc_auc))

    clf = LogisticRegression(random_state=49, solver='lbfgs', multi_class='multinomial')
    clf.fit(myHumanDataReplaced, myLabels)
    mousePred = clf.predict_proba(myMouseDataReplaced)
    fpr, tpr, thresholds = roc_curve(myMouseLabels, mousePred[:,1])
    print(clf.coef_)
    roc_auc = auc(fpr, tpr)
    panel2.plot(fpr, tpr, color='r', label='Logistic Regression (AUC = %0.3f)' % (roc_auc))

    clf = SVC(probability=True, gamma='auto', kernel='linear')
    clf.fit(myHumanDataReplaced, myLabels)
    mousePred = clf.predict_proba(myMouseDataReplaced)
    fpr, tpr, thresholds = roc_curve(myMouseLabels, mousePred[:,1])
    roc_auc = auc(fpr, tpr)
    panel2.plot(fpr, tpr, color='y', label='Support Vector Machine (AUC = %0.3f)' % (roc_auc))

    panel2.set_xlabel("False Positive Rate", fontsize=18)
    panel2.set_ylabel("True Positive Rate", fontsize=18)

    panel1.set_xlim([-0.01,1.01])
    panel1.set_ylim([-0.01,1.01])
    panel2.set_xlim([-0.01,1.01])
    panel2.set_ylim([-0.01,1.01])

    panel1.text(0.01, 1.03, "A", ha='center', va='bottom', fontsize=32)
    panel2.text(0.01, 1.03, "B", ha='center', va='bottom', fontsize=32)

    panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=20)
    panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=20)
    panel1.legend(loc="lower right", fontsize=16)
    panel2.legend(loc="lower right", fontsize=16)


    activeDict, inactiveDict = plotMouseNew()

    panel3 = plt.axes([panel_width*2+left_right*3, 0.57, panel_width/1.1, panel_height/2.47], frameon=True)
    panel4 = plt.axes([panel_width*2+left_right*3, above_below, panel_width/1.1, panel_height/2.47], frameon=True)

    panel3.hist(activeDict,edgecolor='black', linewidth=1, bins=np.arange(11)-0.5)
    panel4.hist(inactiveDict,edgecolor='black', linewidth=1, bins=np.arange(11)-0.5)
    panel3.set_xlabel("Actively Transcribed Tissues", fontsize=18)
    panel3.set_ylabel("Genes Predicted Active", fontsize=18)
    panel4.set_xlabel("Actively Transcribed Tissues", fontsize=18)
    panel4.set_ylabel("Genes Predicted Inactive", fontsize=18)
    panel3.set_ylim([0,103])
    panel4.set_ylim([0,103])
    panel3.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=20)
    panel4.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=20)
    panel3.text(-1, 103*1.03, 'C', ha='center', va='bottom', fontsize=32)
    panel4.text(-1, 103*1.03, 'D', ha='center', va='bottom', fontsize=32)


    plt.savefig('test_combined.pdf', dpi=700)
    plt.close()


def histToList(myList):
    myReturn = [0]*(max(myList)+1)
    for i in myList:
        myReturn[i] += 1
    return(myReturn)


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

    NCNameToPred = {}
    for line in open('Mmusc10tRNAClassificationsNewFixedNoSegDups.txt'):
        splitLine = (line.strip()).split()
        NCNameToPred[splitLine[0]] = float(splitLine[1])

    activity = []
    pred = []
    activeDict = []
    inactiveDict = []

    NCToAboveZero = {}
    for line in open('mouse_liver_1.bed'):
        splitLine = (line.strip()).split('\t')
        if not splitLine[0] == 'tRNA':
            if int(splitLine[3]) > 0:
                NCToAboveZero[splitLine[0]] = True
            else:
                NCToAboveZero[splitLine[0]] = False

    activeZero = 0
    predActiveZero = 0
    for tRNA in tRNAToNCName:
        if tRNAToNCName[tRNA] in NCNameToPred:
            activity.append(tRNAToNumActive[tRNA])
            pred.append(NCNameToPred[tRNAToNCName[tRNA]])
            if NCNameToPred[tRNAToNCName[tRNA]] < 0.0:
                inactiveDict.append(tRNAToNumActive[tRNA])
            elif NCNameToPred[tRNAToNCName[tRNA]] > 0.0:
                activeDict.append(tRNAToNumActive[tRNA])

            if tRNAToNCName[tRNA] in NCToAboveZero and NCToAboveZero[tRNAToNCName[tRNA]] == False:
                if tRNAToNumActive[tRNA] > 0:
                    activeZero += 1
                    if NCNameToPred[tRNAToNCName[tRNA]] > 0.0:
                        predActiveZero += 1

    print(activeZero)
    print(predActiveZero)
    print(activeDict.count(1))
    print(len(activeDict)+len(inactiveDict))


    # print(inactiveDict)
    # print(activeDict)

    activeList = histToList(activeDict)
    inactiveList = histToList(inactiveDict)
    inactiveList.append(0)
    activeRange = list(range(0,len(activeList)))
    inactiveRange = []
    for k in activeRange:
        inactiveRange.append(k+0.5)

    print(activeDict)
    print(inactiveDict)
    print(activeRange)
    print(inactiveRange)
    return(activeDict, inactiveDict)





def reorder(myList, myOrder):
    myReturn = []
    for k in myOrder:
        myReturn.append(myList[k])
    return(myReturn)

def makeFloat(myList):
    myReturn = []
    for k in myList:
        if not k == '?' and not 'tRNA' in k:
            myReturn.append(float(k))
        elif k == '?':
            myReturn.append(np.nan)
        else:
            myReturn.append(k)
    return(myReturn)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\n'.join(newList)

def getScore(list1,list2):
    myTotal = 0.0
    for i in range(0,len(list1)):
        if ((list1[i])[0] > (list1[i])[1] and int(list2[i]) == 0) or ((list1[i])[1] > (list1[i])[0] and int(list2[i]) == 1):
            myTotal += 1.0
    return(myTotal/float(len(list2)))

def confusionMatrix(pred,real):
    myReturn = [[0.0, 0.0],[0.0, 0.0]]
    for i in range(0,len(pred)):
        if int(pred[i]) == 0 and int(real[i]) == 0:
            (myReturn[0])[0] += 1
        elif int(pred[i]) == 0 and int(real[i]) == 1:
            (myReturn[0])[1] += 1
        elif int(pred[i]) == 1 and int(real[i]) == 0:
            (myReturn[1])[0] += 1
        elif int(pred[i]) == 1 and int(real[i]) == 1:
            (myReturn[1])[1] += 1
    return(myReturn)

def getAccuracy(cM):
    return( (((cM[0])[0]+(cM[1])[1]) / ((cM[0])[0]+(cM[0])[1]+(cM[1])[0]+(cM[1])[1])) * 100.0 )

def getFirst(myList):
    myReturn = []
    for k in myList:
        myReturn.append(k[0])
    return(myReturn)



def main():
    buildClassifier()
    
if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit