#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 14:47:47 2018

@author: neeks
"""

import numpy as np
from scipy.io import wavfile

path = '/Users/neeks/data_temp/'
#/home/neeks/mnt_1/work/neeks/spkrChange/dBase/libriSpeech/train-clean-100-sel/'

#fID = [0374, 2843, 5456, 7447, 7505] # this is a list

fID = ['0374', '2843', '5456', '7447', '7505'] # this is a list

for i in range(len(fID)):
    fname = path + str(fID[i]) + '/' + str(fID[i]) + '_trans.txt'
    print('Reading file:' + fname)
    
    with open(fname) as f:
        content = f.readlines()
   
    # get the lines
    content = [x.strip() for x in content]
    
    # get the words in each line
    word_rate = []
    for j in range(len(content)):
        words = content[j].split(' ')
        nos = len(words)-1

        # read wave file to find duration
        wav_fname = path + str(fID[i]) + '/' + words[0] + '.wav'
        sampFreq, snd = wavfile.read(wav_fname)
        snd = snd/np.max(snd)
        dur = float(len(snd))/sampFreq
#         print(dur)
        word_rate.append(float(nos)/float(dur))
    print('ID: '+str(fID[i])+' mean word rate = '+str(np.mean(word_rate)) +' words-per-s')
    print('with std. dev = '+str(np.std(word_rate)) +' words-per-s')