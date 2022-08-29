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
allSubjects = studyparams.EPHYS_MICE
#allSubjects = studyparams.TEST_MOUSE



for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse
    dbPath = os.path.join(databaseDir, f'{subject}_paspeech_speech_pval.h5')

    celldb = celldatabase.load_hdf(dbPath)

    nCells = len(celldb)

    newdbPath = os.path.join(databaseDir , f'{subject}_paspeech_speech_tuning.h5')
    celldb['ftSelectivityVotMaxPvalOnset'] = np.nan
    celldb['ftSelectivityVotMaxPvalSustain'] = np.nan
    celldb['ftSelectivityVotMinPvalOnset'] = np.nan
    celldb['ftSelectivityVotMinPvalSustain'] = np.nan

    celldb['votSelectivityFtMaxPvalOnset'] = np.nan
    celldb['votSelectivityFtMaxPvalSustain'] = np.nan
    celldb['votSelectivityFtMinPvalOnset'] = np.nan
    celldb['votSelectivityFtMinPvalSustain'] = np.nan

    periodsName = ['respOnset', 'respSustained']
    allPeriods = [ [0, 0.12] , [0.12, 0.24] ]
    periodDuration = [x[1]-x[0] for x in allPeriods]
    N_FT = 4 # HARDCODED
    N_VOT = 4 #HARDCODED
    N_SPEECH = 12 #HARDCODED


    correctedAlpha = 0.05/N_SPEECH
    celldbResp = celldb[(celldb.speechMinPvalOnset < correctedAlpha) |
                        (celldb.speechMinPvalSustain < correctedAlpha)]


    indCell = -1
    for indRow, dbRow in celldbResp.iterrows():
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


        FTParamsEachTrial = bdata['targetFTpercent']
        possibleFTParams = np.unique(FTParamsEachTrial)
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        nFT = len(possibleFTParams)
        nVOT = len(possibleVOTParams)
        minVOT = [np.min(possibleVOTParams)]
        maxVOT = np.max(possibleVOTParams)
        minFT = np.min(possibleFTParams)
        maxFT = np.max(possibleFTParams)

        trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)

        FTSelectivity_VOTMin_Pval = np.full(len(allPeriods), np.nan)
        FTSelectivity_VOTMax_Pval = np.full(len(allPeriods), np.nan)
        VOTSelectivity_FTMin_Pval = np.full(len(allPeriods), np.nan)
        VOTSelectivity_FTMax_Pval = np.full(len(allPeriods), np.nan)

        # calculate selectvitiy to FT
        for indPeriod, period in enumerate(allPeriods):
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial,
                                                                     period)
            nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

            # -- Calculate nSpikes for each FT @ VOTmin and VOTmax, and each VOT @ FTmin and FTmax

            nSpikesEachFT_VOTmin = []
            nSpikesEachFT_VOTmax = []
            nSpikesEachVOT_FTmin = []
            nSpikesEachVOT_FTmax = []

            for indcond, thisCond in enumerate(possibleFTParams):
                nSpikesEachFT_VOTmin.append(nSpikesEachTrial[trialsEachCond[:,indcond,0]])
                nSpikesEachFT_VOTmax.append(nSpikesEachTrial[trialsEachCond[:,indcond,3]])
                nSpikesEachVOT_FTmin.append(nSpikesEachTrial[trialsEachCond[:,0,indcond]])
                nSpikesEachVOT_FTmax.append(nSpikesEachTrial[trialsEachCond[:,3,indcond]])

            if np.all(spikeCountMat == 0):
                kStat = None
                pValKruskal = 1
                print(f'{oneCell} no spikes for FT-VOT session')
            else:

                try:
                    kStat, pValKruskalFT_VOTmin = stats.kruskal(*nSpikesEachFT_VOTmin)
                except ValueError:
                    pValKruskalFT_VOTmin = 1
                FTSelectivity_VOTMin_Pval[indPeriod] = pValKruskalFT_VOTmin

                try:
                    kStat, pValKruskalFT_VOTmax = stats.kruskal(*nSpikesEachFT_VOTmax)
                except ValueError:
                    pValKruskalFT_VOTmax = 1
                FTSelectivity_VOTMax_Pval[indPeriod] = pValKruskalFT_VOTmax

                try:
                    kStat, pValKruskalVOT_FTmin = stats.kruskal(*nSpikesEachVOT_FTmin)
                except ValueError:
                    pValKruskalVOT_FTmin = 1
                VOTSelectivity_FTMin_Pval[indPeriod] = pValKruskalVOT_FTmin

                try:
                    kStat, pValKruskalVOT_FTmax = stats.kruskal(*nSpikesEachVOT_FTmax)
                except ValueError:
                    pValKruskalVOT_FTmax = 1
                VOTSelectivity_FTMax_Pval[indPeriod] = pValKruskalVOT_FTmax


                '''
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
                '''

    celldb.at[indRow, 'ftSelectivityVotMaxPvalOnset'] = FTSelectivity_VOTMax_Pval[0]
    celldb.at[indRow, 'ftSelectivityVotMaxPvalSustain'] = FTSelectivity_VOTMax_Pval[1]
    celldb.at[indRow, 'ftSelectivityVotMinPvalOnset'] = FTSelectivity_VOTMin_Pval[0]
    celldb.at[indRow, 'ftSelectivityVotMinPvalSustain'] = FTSelectivity_VOTMin_Pval[1]

    celldb.at[indRow, 'votSelectivityFtMaxPvalOnset'] = VOTSelectivity_FTMax_Pval[0]
    celldb.at[indRow, 'votSelectivityFtMaxPvalSustain'] = VOTSelectivity_FTMax_Pval[1]
    celldb.at[indRow, 'votSelectivityFtMinPvalOnset'] = VOTSelectivity_FTMin_Pval[0]
    celldb.at[indRow, 'votSelectivityFtMinPvalSustain'] = VOTSelectivity_FTMin_Pval[1]


    celldatabase.save_hdf(celldb, newdbPath)
