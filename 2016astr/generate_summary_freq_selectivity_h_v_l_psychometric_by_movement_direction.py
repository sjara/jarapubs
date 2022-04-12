'''
Quantify high vs low frequency selectivity fixing one movement direction at a time, using all good cells recorded in the 2afc psychometric task. For all the non-duplicated good cells that are in the striatum, look at their response each trial each frequency for trials either with right or left subsequent choice.
Lan Guo 20180125  
'''

import os
import sys
import numpy as np
import pandas as pd
import importlib
from jaratoolbox import settings
import figparams
from jaratoolbox import loadbehavior
from jaratoolbox import loadopenephys
from jaratoolbox import spikesanalysis
from jaratoolbox import behavioranalysis
import scipy.stats as stats
import pdb

FIGNAME = 'sound_freq_selectivity'
outputDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
dataDir = outputDir

paradigm = '2afc'
scriptFullPath = os.path.realpath(__file__)

qualityList = [1,6]
ISIcutoff = 0.02
numOfFreqs = 6
maxNumOfTrials = 300 # a big number to make sure allocate enough space to store firing rate for all trials of each frequency
EPHYS_SAMPLING_RATE = 30000.0
soundTriggerChannel = 0
# -- Access mounted behavior and ephys drives for psycurve and switching mice -- #
BEHAVIOR_PATH = settings.BEHAVIOR_PATH_REMOTE
EPHYS_PATH = settings.EPHYS_PATH_REMOTE

if not os.path.ismount(BEHAVIOR_PATH):
    os.system('sshfs -o idmap=user jarauser@jarahub:/data/behavior/ {}'.format(BEHAVIOR_PATH))

if not os.path.ismount(EPHYS_PATH):
    os.system('sshfs -o idmap=user jarauser@jarastore:/data2016/ephys/ {}'.format(EPHYS_PATH))


# -- Read in databases storing all measurements from psycurve mice -- #
psychometricFilePath = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME)
psychometricFileName = 'all_cells_all_measures_waveform_psychometric.h5'
psychometricFullPath = os.path.join(psychometricFilePath,psychometricFileName)
allcells_psychometric = pd.read_hdf(psychometricFullPath,key='psychometric')

goodcells_psychometric = (allcells_psychometric.cellQuality.isin(qualityList)) & (allcells_psychometric.ISI <= ISIcutoff)
cellInStr =  (allcells_psychometric.cellInStr==1)
keepAfterDupTest = allcells_psychometric.keep_after_dup_test
cellSelector = goodcells_psychometric & cellInStr & keepAfterDupTest  #Boolean array
cellsToPlot = allcells_psychometric[cellSelector]

responseFilename = 'response_each_freq_each_cell_psycurve_2afc.npz'
responseFullPath = os.path.join(dataDir,responseFilename)
responseEachCellEachFreq = np.load(responseFullPath)

cellsToPlot = cellsToPlot.reset_index()

responseEachFreqEachCellLeftChoice = np.ma.masked_array(np.empty((len(cellsToPlot),maxNumOfTrials,numOfFreqs)),mask=np.zeros((len(cellsToPlot),maxNumOfTrials,numOfFreqs)))
responseEachFreqEachCellRightChoice = np.ma.masked_array(np.empty((len(cellsToPlot),maxNumOfTrials,numOfFreqs)),mask=np.zeros((len(cellsToPlot),maxNumOfTrials,numOfFreqs)))

for indc,cell in cellsToPlot.iterrows(): 
    print 'retrieving data for cell', indc
    animalName = cell['animalName']
    allcellsFileName = 'allcells_'+animalName+'_quality' #This is specific to Billy's final allcells files after adding cluster quality info 
    sys.path.append(settings.ALLCELLS_PATH)
    allcells = importlib.import_module(allcellsFileName)

    behavSession = cell['behavSession']
    tetrode = cell['tetrode']
    cluster = cell['cluster']
    ### Using cellDB methode to find the index of this cell in the cellDB ###
    cellIndex = allcells.cellDB.findcell(animalName, behavSession, tetrode, cluster)
    oneCell = allcells.cellDB[cellIndex]
    
    ephysSession = oneCell.ephysSession
    ## Get behavior data associated with 2afc session ###
    behavFileName = '{0}_{1}_{2}.h5'.format(animalName,paradigm,behavSession)
    behavFile = os.path.join(BEHAVIOR_PATH,animalName,behavFileName)
    bdata = loadbehavior.BehaviorData(behavFile,readmode='full')

    ### Get events data ###
    fullEventFilename=os.path.join(EPHYS_PATH, animalName, ephysSession, 'all_channels.events')
    eventData = loadopenephys.Events(fullEventFilename)
    ##### Get event onset times #####
    eventData.timestamps = np.array(eventData.timestamps)/EPHYS_SAMPLING_RATE #hard-coded ephys sampling rate!!

    eventOnsetTimes=np.array(eventData.timestamps)
    soundOnsetEvents = (eventData.eventID==1) & (eventData.eventChannel==soundTriggerChannel)
    soundOnsetTimes = eventOnsetTimes[soundOnsetEvents]
    soundOnsetTimeBehav = bdata['timeTarget']
    # -- Check to see if ephys and behav recordings have same number of trials, remove missing trials from behav file -- #
    # Find missing trials
    missingTrials = behavioranalysis.find_missing_trials(soundOnsetTimes,soundOnsetTimeBehav)
    # Remove missing trials
    bdata.remove_trials(missingTrials)
    possibleFreq = np.unique(bdata['targetFrequency'])
    numFreqs = len(possibleFreq)

    responseEachFreqLeftChoice = [] #store sound-evoked response (0-100ms) for each trial with left choice
    responseEachFreqRightChoice = [] #store sound-evoked response (0-100ms) for each trial with right choice
    for indf, freq in enumerate(possibleFreq):
        # -- Only use valid trials of one frequency to estimate response index -- #
        rightward = bdata['choice']==bdata.labels['choice']['right']
        leftward = bdata['choice']==bdata.labels['choice']['left'] 
        oneFreqTrials = (bdata['targetFrequency'] == freq) & bdata['valid'].astype('bool')
        oneFreqTrialsLeftChoice = leftward[oneFreqTrials] 
        oneFreqTrialsRightChoice = rightward[oneFreqTrials] 
        # need to figure out which trials of the oneFreqTrials are left, which are right
        oneFreqResps = responseEachCellEachFreq[indc,:,indf].compressed() 
        oneFreqRespsLeftChoice = oneFreqResps[oneFreqTrialsLeftChoice]
        oneFreqRespsRightChoice = oneFreqResps[oneFreqTrialsRightChoice]
        #responseEachFreqLeftChoice.append(oneFreqRespsLeftChoice)
        #responseEachFreqRightChoice.append(oneFreqRespsRightChoice)
        numOfTrialsLeft = len(oneFreqRespsLeftChoice)
        responseEachFreqEachCellLeftChoice[indc,:numOfTrialsLeft,indf] = oneFreqRespsLeftChoice
        responseEachFreqEachCellLeftChoice.mask[indc,numOfTrialsLeft:,indf] = True #Mask the extra trials not occupied by data
        numOfTrialsRight = len(oneFreqRespsRightChoice)
        responseEachFreqEachCellRightChoice[indc,:numOfTrialsRight,indf] = oneFreqRespsRightChoice
        responseEachFreqEachCellRightChoice.mask[indc,numOfTrialsRight:,indf] = True #Mask the extra trials not occupied by data

# -- Save psth intermediate data -- #
if not os.path.exists(outputDir):
    os.mkdir(outputDir)


# Have to save the two masked arrays individually because:
# NotImplementedError: MaskedArray.tofile() not implemented yet.
# Cannot directly save masked arrays in npz
# responseEachFreqEachCell=responseEachFreqEachCell, baselineEachFreqEachCell=baselineEachFreqEachCell,
responseEachFreqEachCellLeftChoiceFile = 'response_each_freq_each_cell_psycurve_left_choice.npz'
responseEachFreqEachCellRightChoiceFile = 'response_each_freq_each_cell_psycurve_right_choice.npz'

responseEachFreqEachCellLeftChoice.dump(os.path.join(outputDir,responseEachFreqEachCellLeftChoiceFile))
responseEachFreqEachCellRightChoice.dump(os.path.join(outputDir,responseEachFreqEachCellRightChoiceFile))
