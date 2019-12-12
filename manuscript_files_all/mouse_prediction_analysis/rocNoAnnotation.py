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
    labelsDict['SecStruct'] = "tRNAscan-SE 2'\nStructure Score"
    labelsDict['HMMScore'] = "tRNAscan-SE\nHMM Score"
    labelsDict['IsoScore'] = "tRNAscan-SE Isotype\nBit Score"
    labelsDict['tRNA10kb'] = 'tRNA Genes within 10 Kilobases'
    labelsDict['Prot75kb'] = 'Exons within 75 Kilobases'
    labelsDict['TTTT'] = 'Distance to Nearest TTTT\nTranscription Termination\nSequence'
    labelsDict['Codon'] = 'tRNAs Corresponding\nto the Same Codon'
    labelsDict['MFE'] = 'Constrained Minimum\nFree Energy'

    myHumanData = []
    myLabels = []
    myHumanNames = []
    imp_mean = SimpleImputer(missing_values=np.nan, strategy='mean')
    for line in open('humanNoAnnotation.tsv'):
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
    for line in open('mouseNoAnnotation.tsv'):
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
    open('humanCVPredictionsChangedNoAnnotation.txt', 'w').write(myOutString)



    print(clf.feature_importances_)

    # cM = confusionMatrix(cvPredictions, np.asarray(myLabels))
    # print(cM)
    # print(getAccuracy(cM))


    fig_width = 14
    fig_height = 7
    plt.figure(figsize=(fig_width, fig_height))
    panel_width = 0.4
    panel_height = 0.8
    panel_total_height = (panel_height*1)
    extra_y_space = 1 - panel_total_height
    above_below = extra_y_space/2
    panel_total_width = (panel_width*2)
    extra_x_space = 1 - panel_total_width
    left_right = extra_y_space/3
    panel1 = plt.axes([left_right, (1-panel_height)/2, panel_width, panel_height], frameon=True)
    panel2 = plt.axes([1-panel_width-left_right, (1-panel_height)/2, panel_width, panel_height], frameon=True)

    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    fpr, tpr, thresholds = roc_curve(myLabels, cvPredictions[:,1])
    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='b', label='Random Forest (AUC = %0.3f)' % (roc_auc))

    clf = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    fpr, tpr, thresholds = roc_curve(myLabels, cvPredictions[:,1])
    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='r', label='Logistic Regression (AUC = %0.3f)' % (roc_auc))
    print(clf.coef_)

    clf = SVC(probability=True, gamma='auto', kernel='linear')
    clf.fit(myHumanDataReplaced, myLabels)
    cvPredictions = cross_val_predict(clf, myHumanDataReplaced, myLabels, cv=10, method='predict_proba')
    fpr, tpr, thresholds = roc_curve(myLabels, cvPredictions[:,1])
    roc_auc = auc(fpr, tpr)
    panel1.plot(fpr, tpr, color='y', label='Support Vector Machine (AUC = %0.3f)' % (roc_auc))

    panel1.set_xlabel("False Positive Rate", fontsize=15)
    panel1.set_ylabel("True Positive Rate", fontsize=15)

    clf = RandomForestClassifier(n_estimators=250, max_depth=4, random_state=49, oob_score=True, n_jobs=8, min_samples_split=2)
    clf.fit(myHumanDataReplaced, myLabels)
    mousePred = clf.predict_proba(myMouseDataReplaced)
    fpr, tpr, thresholds = roc_curve(myMouseLabels, mousePred[:,1])
    roc_auc = auc(fpr, tpr)
    panel2.plot(fpr, tpr, color='b', label='Random Forest (AUC = %0.3f)' % (roc_auc))

    clf = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial')
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

    panel2.set_xlabel("False Positive Rate", fontsize=15)
    panel2.set_ylabel("True Positive Rate", fontsize=15)

    panel1.set_xlim([-0.01,1.01])
    panel1.set_ylim([-0.01,1.01])
    panel2.set_xlim([-0.01,1.01])
    panel2.set_ylim([-0.01,1.01])

    panel1.text(0.01, 1.03, "A", ha='center', va='bottom', fontsize=32)
    panel2.text(0.01, 1.03, "B", ha='center', va='bottom', fontsize=32)

    panel1.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=15)
    panel2.tick_params(bottom='on', labelbottom='on',\
               left='on', labelleft='on', \
               right='off', labelright='off',\
               top='off', labeltop='off', labelsize=15)
    panel1.legend(loc="lower right", fontsize=15)
    panel2.legend(loc="lower right", fontsize=15)
    plt.savefig('cvROCChangedNoAnnotation.png', dpi=700)
    plt.close()





    # print(cross_validate(clf, myHumanDataReplaced, myLabels, cv=10, return_train_score=True, return_estimator=True)['estimator'])

    # a_train, a_test, b_train, b_test = train_test_split(myHumanDataReplaced, myLabels, test_size=0.2, random_state=42)
    
    # clf.fit(a_train, b_train)
    # print(getScore(clf.predict_proba(a_test),np.asarray(b_test)))


    # clf = RandomForestClassifier(n_estimators=1000, max_depth=5, random_state=3, oob_score=True, n_jobs=8, min_samples_split=2)
    # clf.fit(myHumanDataReplaced, myLabels)
    # print(clf.score(myMouseDataReplaced,myMouseLabels))
    # cM = confusionMatrix(clf.predict(myMouseDataReplaced),myMouseLabels)
    # print(cM)
    # print(getAccuracy(cM))
    # print(getScore(clf.predict_proba(myMouseDataReplaced),myMouseLabels))

    # clf = LogisticRegression(random_state=0, solver='lbfgs', multi_class='multinomial').fit(myHumanData / np.std(myHumanData, 0), myLabels)
    # print(clf.coef_)
    # print(clf.score(myMouseData,myMouseLabels))

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