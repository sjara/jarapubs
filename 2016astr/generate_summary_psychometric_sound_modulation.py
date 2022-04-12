'''
Generate and store intermediate data for plot showing response to same stimulus being modulated by different choice of action in the psychometric curve tasks. 
Output npz file contains animalName, modulated(boolean vector indicating whether each cell was significantly modulated), and modulationIndex(modulation index for each cell).
Lan Guo20170103
'''
# TO DO
# Include in the database information about what script creates it
# Save cell databases on jarahub

import os
import numpy as np
import pandas as pd
from jaratoolbox import settings
import figparams

scriptFullPath = os.path.realpath(__file__)
qualityList = [1,6]
maxZThreshold = 3

# -- Read in databases storing all measurements from switching mice -- #
psychometricFilePath = os.path.join(settings.FIGURESDATA,figparams.STUDY_NAME)
psychometricFileName = 'all_cells_all_measures_waveform_psychometric.h5'
psychometricFullPath = os.path.join(psychometricFilePath,psychometricFileName)
allcells_psychometric = pd.read_hdf(psychometricFullPath,key='psychometric')

# -- Get intermediate data relevant to this subfigure of sound modulation -- #
## Plot all good cells from switcing task ##
goodcells = (allcells_psychometric.cellQuality.isin(qualityList)) & (allcells_psychometric.ISI <= 0.02)
allcells = allcells_psychometric[goodcells & (allcells_psychometric.cellInStr==1)] #just look at the good cells THAT ARE IN STR 
strongerSoundResMid1 = abs(allcells.maxZSoundMid1) > abs(allcells.maxZSoundMid2)
modIStrongerSoundRes = np.r_[allcells.modIndexMid1[strongerSoundResMid1].values,allcells.modIndexMid2[~strongerSoundResMid1].values]
modSigStrongerSoundRes = np.r_[allcells.modSigMid1[strongerSoundResMid1].values,allcells.modSigMid2[~strongerSoundResMid1].values] 
#sigMod = np.array(((allcells.modSigMid1<=0.05)|(allcells.modSigMid2<=0.05)), dtype=bool)
sigMod = np.array((modSigStrongerSoundRes <= 0.05), dtype=bool)
dataToPlot = {'modulated':sigMod,'modulationIndex':modIStrongerSoundRes}

### Save data ###
outputDir = os.path.join(settings.FIGURESDATA, figparams.STUDY_NAME)
outputFile = 'summary_psychometric_sound_modulation_all_good_cells.npz'
outputFullPath = os.path.join(outputDir,outputFile)
np.savez(outputFullPath, sourcePsychometric=psychometricFilePath, goodCellQuality=qualityList, script=scriptFullPath, **dataToPlot)


## Plot only the good cells responsive to middle freq from switcing task ##
responsiveMidFreqs = (abs(allcells_psychometric.maxZSoundMid1)>=maxZThreshold) | (abs(allcells_psychometric.maxZSoundMid2)>=maxZThreshold)
allcellsResponsive = allcells_psychometric[goodcells & responsiveMidFreqs & (allcells_psychometric.cellInStr==1)] # good quality cells IN STR & sound responsive to middle frequency
strongerSoundResMid1 = abs(allcellsResponsive.maxZSoundMid1) > abs(allcellsResponsive.maxZSoundMid2)
modIStrongerSoundRes = np.r_[allcellsResponsive.modIndexMid1[strongerSoundResMid1].values,allcellsResponsive.modIndexMid2[~strongerSoundResMid1].values]
modSigStrongerSoundRes = np.r_[allcellsResponsive.modSigMid1[strongerSoundResMid1].values,allcellsResponsive.modSigMid2[~strongerSoundResMid1].values] 
#sigMod = np.array(((allcellsResponsive.modSigMid1<=0.05)|(allcellsResponsive.modSigMid2<=0.05)), dtype=bool)
sigMod = np.array((modSigStrongerSoundRes <= 0.05), dtype=bool)
dataToPlot = {'modulated':sigMod,'modulationIndex':modIStrongerSoundRes}

### Save data ###
outputDir = os.path.join(settings.FIGURESDATA, figparams.STUDY_NAME)
outputFile = 'summary_psychometric_sound_modulation_good_cells_responsive_midfreq.npz'
outputFullPath = os.path.join(outputDir,outputFile)
np.savez(outputFullPath, sourcePsychometric=psychometricFilePath, maxZThreshold=maxZThreshold, goodCellQuality=qualityList, script=scriptFullPath, **dataToPlot)
