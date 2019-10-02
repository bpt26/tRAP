#!/usr/bin/env python3
# Name: Bryan Thornlow
# Date: 4/10/2018
# findSegDups.py

import sys
import os
import time
import random
import numpy
import gzip
import math

##########################
##### MAIN FUNCTIONS #####
##########################

def removeSegDupsFromAlignment():
	abbs = ['Ananc2','Btaur8','Casia1']
	abbs += ['Cfami3','Chirc1','Clani1']
	abbs += ['Cporc3','Dnove3','Ecaba2']
	abbs += ['Eeuro2','Efusc1','Ggori5']
	abbs += ['Hglab2','Hsapi38','Jjacu1']
	abbs += ['Mmarm2','Mmula8','Mmuri3']
	abbs += ['Mmusc10','Mochr1','Mputo1']
	abbs += ['Ocuni2','Odegu1','Oorca1']
	abbs += ['Palec1','Ppygm3','Ptrog5']
	abbs += ['Rnorv6','Sscro11']

	toRemove = []
	for abb in abbs:
		for line in open('segDupFiles/'+abb.upper()+'SegDups.txt'):
			toRemove.append(abb+'-'+line.strip())
	print(toRemove)

	toKeep = joiner(sorted(abbs))+'\n'
	for line in open('oneToOneMapAllSpeciesAugmentMacaque.txt'):
		splitLine = (line.strip()).split('\t')
		if not splitLine[0] == 'Ananc2':
			lineKeep = True
			if splitLine.count('N/A') == 28:
				for k in splitLine:
					if k != 'N/A':
						if k in toRemove:
							lineKeep = False
			if lineKeep == True:
				toKeep += joiner(splitLine)+'\n'

	open('FINAL_ORTHOLOG_SET.txt', 'w').write(toKeep)

def joiner(entry):
    newList = []
    for k in entry:
        newList.append(str(k))
    return '\t'.join(newList)

def main():
    removeSegDupsFromAlignment()


if __name__ == "__main__":
    """
    Calls main when program is run by user.
    """
    main();
    raise SystemExit



