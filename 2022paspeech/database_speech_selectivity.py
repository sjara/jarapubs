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
#allSubjects = studyparams.TEST_MOUSE
allSubjects = ['feat010']



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


    celldb['selectivityIndexFT_VOTmin'] = np.nan
    celldb['selectivityIndexFT_VOTmax'] = np.nan
    celldb['selectivityIndexVOT_FTmin'] = np.nan
    celldb['selectivityIndexVOT_FTmin'] = np.nan

    celldb['selectivityIndex2FT_VOTmin'] = np.nan
    celldb['selectivityIndex2FT_VOTmax'] = np.nan
    celldb['selectivityIndex2VOT_FTmin'] = np.nan
    celldb['selectivityIndex2VOT_FTmin'] = np.nan


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

        FTParamsEachTrial = bdata['targetFTpercent']
        possibleFTParams = np.unique(FTParamsEachTrial)
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        nFT = len(possibleFTParams)
        nVOT = len(possibleVOTParams)

        if len(FTParamsEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(FTParamsEachTrial)]
            print(f'[{indRow}] Warning! BehavTrials ({len(freqEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')

        trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)

        FTSelectivity_VOTMin_Pval = np.full(len(allPeriods), np.nan)
        FTSelectivity_VOTMax_Pval = np.full(len(allPeriods), np.nan)
        VOTSelectivity_FTMin_Pval = np.full(len(allPeriods), np.nan)
        VOTSelectivity_FTMax_Pval = np.full(len(allPeriods), np.nan)

        avgSpikesEachFT_VOTmin = np.full([4, len(allPeriods)], np.nan)
        avgSpikesEachFT_VOTmax = np.full([4, len(allPeriods)], np.nan)
        avgSpikesEachVOT_FTmin = np.full([4, len(allPeriods)], np.nan)
        avgSpikesEachVOT_FTmax = np.full([4, len(allPeriods)], np.nan)

        selectivityIndexFT_VOTmin = np.zeros([len(allPeriods)])
        selectivityIndexFT_VOTmax = np.zeros([len(allPeriods)])
        selectivityIndexVOT_FTmin = np.zeros([len(allPeriods)])
        selectivityIndexVOT_FTmax = np.zeros([len(allPeriods)])

        selectivityIndex2FT_VOTmin = np.zeros([len(allPeriods)])
        selectivityIndex2FT_VOTmax = np.zeros([len(allPeriods)])
        selectivityIndex2VOT_FTmin = np.zeros([len(allPeriods)])
        selectivityIndex2VOT_FTmax = np.zeros([len(allPeriods)])


        # calculate spike times
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        for indPeriod, period in enumerate(allPeriods):

            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial, period)
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


            ## Calculate selectivity index
            for indCond, thisCond in enumerate(possibleFTParams):
                avgSpikesEachFT_VOTmin[indCond, indPeriod] = np.mean(nSpikesEachFT_VOTmin[indCond])
                avgSpikesEachFT_VOTmax[indCond, indPeriod] = np.mean(nSpikesEachFT_VOTmax[indCond])
                avgSpikesEachVOT_FTmin[indCond, indPeriod] = np.mean(nSpikesEachVOT_FTmin[indCond])
                avgSpikesEachVOT_FTmax[indCond, indPeriod] = np.mean(nSpikesEachVOT_FTmax[indCond])
            # Selectivity index 1: (max - min)/(max + min)

            if np.all(spikeCountMat == 0):
                selectivityIndexFT_VOTmin[indPeriod] = 0
                selectivityIndexFT_VOTmax[indPeriod] = 0
                selectivityIndexVOT_FTmin[indPeriod] = 0
                selectivityIndexVOT_FTmax[indPeriod] = 0

                selectivityIndex2FT_VOTmin[indPeriod] = 0
                selectivityIndex2FT_VOTmax[indPeriod] = 0
                selectivityIndex2VOT_FTmin[indPeriod] = 0
                selectivityIndex2VOT_FTmax[indPeriod] = 0

                print(f'{oneCell} no spikes for FT-VOT session')
            else:
                selectivityIndexFT_VOTmin[indPeriod] = (np.max(avgSpikesEachFT_VOTmin[:,indPeriod]) - np.min(avgSpikesEachFT_VOTmin[:,indPeriod])) /(np.max(avgSpikesEachFT_VOTmin[:,indPeriod]) + np.min(avgSpikesEachFT_VOTmin[:,indPeriod]))
                selectivityIndexFT_VOTmax[indPeriod] = (np.max(avgSpikesEachFT_VOTmax[:,indPeriod]) - np.min(avgSpikesEachFT_VOTmax[:,indPeriod])) /(np.max(avgSpikesEachFT_VOTmax[:,indPeriod]) + np.min(avgSpikesEachFT_VOTmax[:,indPeriod]))
                selectivityIndexVOT_FTmin[indPeriod] = (np.max(avgSpikesEachVOT_FTmin[:,indPeriod]) - np.min(avgSpikesEachVOT_FTmin[:,indPeriod]))/(np.max(avgSpikesEachVOT_FTmin[:,indPeriod]) + np.min(avgSpikesEachVOT_FTmin[:,indPeriod]))
                selectivityIndexVOT_FTmax[indPeriod] = (np.max(avgSpikesEachVOT_FTmax[:,indPeriod]) - np.min(avgSpikesEachVOT_FTmax[:,indPeriod]))/(np.max(avgSpikesEachVOT_FTmax[:,indPeriod]) + np.min(avgSpikesEachVOT_FTmax[:,indPeriod]))
                # Selectivity index 2 = (max-min)/min
                selectivityIndex2FT_VOTmin[indPeriod] = (np.max(avgSpikesEachFT_VOTmin[:,indPeriod]) - np.min(avgSpikesEachFT_VOTmin[:,indPeriod]))/ np.min(avgSpikesEachFT_VOTmin[:,indPeriod])
                selectivityIndex2FT_VOTmax[indPeriod] = (np.max(avgSpikesEachFT_VOTmax[:,indPeriod]) - np.min(avgSpikesEachFT_VOTmax[:,indPeriod]))/ np.min(avgSpikesEachFT_VOTmax[:,indPeriod])
                selectivityIndex2VOT_FTmin[indPeriod] = (np.max(avgSpikesEachVOT_FTmin[:,indPeriod]) - np.min(avgSpikesEachVOT_FTmin[:,indPeriod]))/ np.min(avgSpikesEachVOT_FTmin[:,indPeriod])
                selectivityIndex2VOT_FTmax[indPeriod] = (np.max(avgSpikesEachVOT_FTmax[:,indPeriod]) - np.min(avgSpikesEachVOT_FTmax[:,indPeriod]))/ np.min(avgSpikesEachVOT_FTmax[:,indPeriod])




        celldb.at[indRow, 'selectivityIndexFT_VOTmin'] = selectivityIndexFT_VOTmin[0]
        celldb.at[indRow, 'selectivityIndexFT_VOTmax'] = selectivityIndexFT_VOTmax[0]
        celldb.at[indRow, 'selectivityIndexVOT_FTmin'] = selectivityIndexVOT_FTmin[0]
        celldb.at[indRow, 'selectivityIndexVOT_FTmax'] = selectivityIndexVOT_FTmax[0]

        celldb.at[indRow, 'selectivityIndex2FT_VOTmin'] = selectivityIndex2FT_VOTmin[0]
        celldb.at[indRow, 'selectivityIndex2FT_VOTmax'] = selectivityIndex2FT_VOTmax[0]
        celldb.at[indRow, 'selectivityIndex2VOT_FTmin'] = selectivityIndex2VOT_FTmin[0]
        celldb.at[indRow, 'selectivityIndex2VOT_FTmax'] = selectivityIndex2VOT_FTmax[0]


        celldb.at[indRow, 'ftSelectivityVotMaxPvalOnset'] = FTSelectivity_VOTMax_Pval[0]
        celldb.at[indRow, 'ftSelectivityVotMaxPvalSustain'] = FTSelectivity_VOTMax_Pval[1]
        celldb.at[indRow, 'ftSelectivityVotMinPvalOnset'] = FTSelectivity_VOTMin_Pval[0]
        celldb.at[indRow, 'ftSelectivityVotMinPvalSustain'] = FTSelectivity_VOTMin_Pval[1]

        celldb.at[indRow, 'votSelectivityFtMaxPvalOnset'] = VOTSelectivity_FTMax_Pval[0]
        celldb.at[indRow, 'votSelectivityFtMaxPvalSustain'] = VOTSelectivity_FTMax_Pval[1]
        celldb.at[indRow, 'votSelectivityFtMinPvalOnset'] = VOTSelectivity_FTMin_Pval[0]
        celldb.at[indRow, 'votSelectivityFtMinPvalSustain'] = VOTSelectivity_FTMin_Pval[1]


    celldatabase.save_hdf(celldb, newdbPath)
