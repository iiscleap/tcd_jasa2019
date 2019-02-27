#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 15 17:33:28 2018

@author: neeks
"""

import glob # for_listing_filenames_from_directory
import os # for_extracting_filename_from_path
import librosa # audio_processing_library
import h5py # reading/writing h5 files
from yaafelib import *
import numpy as np
import matplotlib.pyplot as plt
# Build a DataFlow object using FeaturePlan

fp = FeaturePlan(sample_rate=16000)
fp.addFeature('mfcc: MFCC blockSize=400 stepSize=160 CepsIgnoreFirstCoeff=1 CepsNbCoeffs=12')
fp.addFeature('mfcc_d1: MFCC blockSize=400 stepSize=160 CepsIgnoreFirstCoeff=1 CepsNbCoeffs=12 > Derivate DOrder=1')
fp.addFeature('mfcc_d2: MFCC blockSize=400 stepSize=160 CepsIgnoreFirstCoeff=1 CepsNbCoeffs=12 > Derivate DOrder=2')

#==============================================================================
# fp.addFeature('engy: Energy blockSize=400  stepSize=160')
# fp.addFeature('envp: EnvelopeShapeStatistics EnDecim=200  blockSize=400  stepSize=160')
# fp.addFeature('lsf: LSF blockSize=400 stepSize=160')
# fp.addFeature('loudness: Loudness FFTLength=0  FFTWindow=Hanning  LMode=Total blockSize=400 stepSize=160')
# fp.addFeature('mel: MelSpectrum blockSize=400  stepSize=160')
# fp.addFeature('psharp: PerceptualSharpness blockSize=400  stepSize=160')
# fp.addFeature('pspread: PerceptualSpread blockSize=400  stepSize=160')
# fp.addFeature('sflat: SpectralFlatness blockSize=400  stepSize=160')
# fp.addFeature('sflux: SpectralFlux blockSize=400  stepSize=160')
# fp.addFeature('sroll: SpectralRolloff blockSize=400  stepSize=160')
# fp.addFeature('sshape: SpectralShapeStatistics blockSize=400  stepSize=160')
# fp.addFeature('sslope: SpectralSlope blockSize=400  stepSize=160')
# fp.addFeature('tshape: TemporalShapeStatistics blockSize=400  stepSize=160')
# fp.addFeature('zcr: ZCR blockSize=400  stepSize=160')
#==============================================================================
engine = Engine()
engine.load(fp.getDataFlow())


# get input metadata
engine.getInputs()
{'audio': {'sampleRate': 16000.0,
           'frameLength': 1,
           'sampleStep': 1,
           'parameters': {'SampleRate': '16000'},
           'size': 1}}

wav_path = './data/libriSpeech/spkrChange/listTest_resp/M_5_spkrs_8_pairs_v0/added_part/'
#==============================================================================
# feats_type = ['mfcc', 'mfcc_d1','mfcc_d2','engy','envp','lsf','loudness','mel','psharp','pspread',
#          'sflat','sflux','sroll','sshape','sslope','tshape','zcr']
# 
#==============================================================================
feats_type = ['mfcc', 'mfcc_d1','mfcc_d2']

store_path = './data/fEATS/libriSpeech/test-listener-wise-data/M_5_spkrs_8_pairs_v0/added_part/'

wav_files = []
for file in glob.glob(wav_path +'LS*.wav'):
    wav_files.append(file)


for i in range(0,len(wav_files)):
    y, sr = librosa.load(wav_files[i],sr=None)
    head, tail = os.path.split(wav_files[i])
    print(str(i+1)+': Extracting features from: ' + tail)
    if sr != 16000:
        y = librosa.resample(y, sr, 16000)
        sr = 16000
        y = y/np.max(abs(y))

    # extract features from a numpy array
    audio = np.reshape(y,(1,len(y)))
    audio = np.float64(audio)
    feats = engine.processAudio(audio)
    
    # save feats
    for j in range(0,len(feats_type)):
        hf = h5py.File(store_path + feats_type[j] + '/' +tail[0:-4] +'.h5','w')
        hf.create_dataset(feats_type[j]+'_framewise',data=feats[feats_type[j]])
        hf.close()
    
        