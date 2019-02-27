#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 13 00:49:24 2018

@author: neeks
"""

import numpy as np
import glob
import os 

annot_path = './data/ibmWatson/output/'
testfile = 'LS_8425_287387_0001_SWT_LS_8425_246962_0030_tChange_6073_ms_ibmWatson_recog.txt'
#file = open(annot_path+testfile, 'r') 
#print file.readline(1)


# create a csv file to store the data
f1 = open("./data/ibmWatson/output/ibmWatson_listest_outputs.csv","w")
f1.write('FNAME,TCHANGE,NSPKRS\n')

# read files into a list
txtfiles = []
for file in glob.glob(annot_path+"*.txt"):
    txtfiles.append(file)

# loop on each file
for i in range(len(txtfiles)):
    # push line by line into an array
    with open(txtfiles[i], "r") as ins:
        fname = os.path.basename(txtfiles[i])
        array = []
        for line in ins:
            array.append(line)  
    # get spkr labels and change instants into array
    index = []
    spkr_label = []
    spkr_start = []
    spkr_end = []
    
    for i in range(len(array)):
        temp = array[i].split()
    #    print(temp)
        if '"speaker":' in temp:
            a = temp[1].strip(',')
            spkr_label.append(int(a))
        if '"from":' in temp:
            a = temp[1].strip(',')
            spkr_start.append(int(1000*float(a)))
        if '"to":' in temp:
            a = temp[1].strip(',')
            spkr_end.append(int(1000*float(a)))
            
    spkr_label = np.asarray(spkr_label)
    spkr_start = np.asarray(spkr_start)
    spkr_end = np.asarray(spkr_end)
            
    x = abs(np.diff(spkr_label))
    indx_change = np.where(x>0)
    nspkr_changes = len(indx_change[0])
    if nspkr_changes>0:
        tchange = spkr_start[indx_change[0][0]+1]
    # write to a csv file
    f1.write("%s,%d,%d \n" % (fname[0:-4],tchange,nspkr_changes))

f1.close()
    

        