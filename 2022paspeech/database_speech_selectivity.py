"""
Calculate speech selectivity
"""

import os
import sys
import studyparams
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
from scipy import signal

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

import studyutils

from importlib import reload
reload(studyutils)


databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#allSubjects = studyparams.EPHYS_MICE
subject = 'feat004'

#for indMouse, thisMouse in enumerate(allSubjects):
#    subject = thisMouse
dbPath = os.path.join(databaseDir, f'{subject}_paspeech_speech_pval.h5')

celldb = celldatabase.load_hdf(dbPath)

nCells = len(celldb)

newdbPath = os.path.join(settings.DATABASE_PATH, f'{subject}_paspeech_speech_tuning.h5')


periodsName = ['respOnset', 'respSustained']
allPeriods = [ [0, 0.12] , [0.12, 0.24] ]
periodDuration = [x[1]-x[0] for x in allPeriods]
N_FT = 4 # HARDCODED
N_VOT = 4 #HARDCODED


correctedAlpha = 0.05/N_FT
celldbRespFT = celldb[(celldb.FTMinPvalOnset < correctedAlpha) |
                    (celldb.FTMinPvalSustain < correctedAlpha)]

celldbRespVOT = celldb[(celldb.VOTMinPvalOnset < correctedAlpha) |
                    (celldb.VOTMinPvalSustain < correctedAlpha)]

celldb['FTSelectivityPvalOnset'] = np.nan
celldb['FTSelectivityPvalSustain'] = np.nan

celldb['VOTSelectivityPvalOnset'] = np.nan
celldb['VOTSelectivityPvalSustain'] = np.nan
#celldb['amMaxSyncRate'] = np.nan

indCell = -1
for indRow, dbRow in celldbRespFT.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('FTVOTBorders')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.4, 0.55]  # In seconds

    #rateEachTrial = bdata['currentFreq']
    # Remove last stim from ephys if not saved in behavior file
    ''' Shouldn't need for speech data.
    if len(rateEachTrial) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(rateEachTrial)]
    '''

    #possibleRate = np.unique(rateEachTrial)
    #nRate = len(possibleRate)
    #trialsEachCond = behavioranalysis.find_trials_each_type(rateEachTrial, possibleRate)

    FTParamsEachTrial = bdata['targetFTpercent']
    possibleFTParams = np.unique(FTParamsEachTrial)
    nFT = len(possibleFTParams)
    trialsEachFTCond = behavioranalysis.find_trials_each_type(FTParamsEachTrial, possibleFTParams)

    VOTParamsEachTrial = bdata['targetVOTpercent']
    possibleVOTParams = np.unique(VOTParamsEachTrial)
    nVOT = len(possibleVOTParams)
    trialsEachVOTCond = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)

    FTSelectivityPval = np.full(len(allPeriods), np.nan)

    # calculate selectvitiy to FT
    for indPeriod, period in enumerate(allPeriods):
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 period)
        nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

        # -- Calculate nSpikes for each rate to test selectivity --
        nSpikesEachFT = []
        for indcond, thisCond in enumerate(possibleFTParams):
            nSpikesEachFT.append(nSpikesEachTrial[trialsEachFTCond[:,indcond]])

        if np.all(spikeCountMat == 0):
            kStat = None
            pValKruskal = None
            print(f'{oneCell} no spikes for FT-VOT session')
        else:
            kStat, pValKruskal = stats.kruskal(*nSpikesEachFT)
            FTSelectivityPval[indPeriod] = pValKruskal

# calculate selectvitiy to VOT
indCell = -1
for indRow, dbRow in celldbRespVOT.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('FTVOTBorders')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.4, 0.55]  # In seconds
    VOTSelectivityPval = np.full(len(allPeriods), np.nan)

    for indPeriod, period in enumerate(allPeriods):
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 period)
        nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

        # -- Calculate nSpikes for each rate to test selectivity --
        nSpikesEachVOT = []
        for indcond, thisCond in enumerate(possibleVOTParams):
            nSpikesEachVOT.append(nSpikesEachTrial[trialsEachVOTCond[:,indcond]])

        if np.all(spikeCountMat == 0):
            kStat = None
            pValKruskal = None
            print(f'{oneCell} no spikes for VOT-VOT session')
        else:
            kStat, pValKruskal = stats.kruskal(*nSpikesEachVOT)
            VOTSelectivityPval[indPeriod] = pValKruskal


celldb.at[indRow, 'FTSelectivityPvalOnset'] = FTSelectivityPval[0]
celldb.at[indRow, 'FTSelectivityPvalSustain'] = FTSelectivityPval[1]
celldb.at[indRow, 'VOTSelectivityPvalOnset'] = VOTSelectivityPval[0]
celldb.at[indRow, 'VOTSelectivityPvalSustain'] = VOTSelectivityPval[1]


#celldb.at[indRow, 'amMaxSyncRate'] = maxSyncRate

#celldatabase.save_hdf(celldb, newdbPath)
