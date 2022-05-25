"""
Calculate AM tuning
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

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse
    dbPath = os.path.join(figuresDataDir, f'{subject}_paspeech_am_pval.h5')

    celldb = celldatabase.load_hdf(dbPath)
    nCells = len(celldb)

    newdbPath = os.path.join('/tmp', f'{subject}_paspeech_am_tuning.h5')


    #responsePeriod = [0, 0.5]
    #respPeriodDuration = responsePeriod[1]-responsePeriod[0]
    periodsName = ['respOnset', 'respSustained']
    allPeriods = [ [0, 0.1] , [0.1, 0.5] ]
    periodDuration = [x[1]-x[0] for x in allPeriods]

    N_RATE = 11 # HARDCODED

    correctedAlpha = 0.05/N_RATE
    celldbResp = celldb[(celldb.amMinPvalOnset < correctedAlpha) |
                        (celldb.amMinPvalSustain < correctedAlpha)]


    #celldbResp = celldbResp.loc[650:]
    #celldbResp = celldbResp.loc[516:]
    #celldbResp = celldbResp.loc[652:]
    #celldbResp = celldbResp.loc[1879:] # Nice sync

    celldb['amSelectivityPvalOnset'] = np.nan
    celldb['amSelectivityPvalSustain'] = np.nan
    celldb['amMaxSyncRate'] = np.nan

    indCell = -1
    for indRow, dbRow in celldbResp.iterrows():
        #dbRow = celldb.loc[570]
        indCell += 1
        oneCell = ephyscore.Cell(dbRow)

        ephysData, bdata = oneCell.load('AM')

        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        timeRange = [-0.4, 0.8]  # In seconds

        rateEachTrial = bdata['currentFreq']
        # Remove last stim from ephys if not saved in behavior file
        if len(rateEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(rateEachTrial)]

        possibleRate = np.unique(rateEachTrial)
        nRate = len(possibleRate)
        trialsEachCond = behavioranalysis.find_trials_each_type(rateEachTrial, possibleRate)

        amSelectivityPvalOnset = np.full(len(allPeriods), np.nan)
        for indPeriod, period in enumerate(allPeriods):
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial,
                                                                     period)
            nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

            # -- Calculate nSpikes for each rate to test selectivity --
            nSpikesEachRate = []
            for indcond, thisCond in enumerate(possibleRate):
                nSpikesEachRate.append(nSpikesEachTrial[trialsEachCond[:,indcond]])
            kStat, pValKruskal = stats.kruskal(*nSpikesEachRate)
            amSelectivityPvalOnset[indPeriod] = pValKruskal

            # -- Evaluate synchronization (during sustained) --
            if indPeriod==1:
                #sustainedSpikes = (spikeTimesFromEventOnset>0.1)&(spikeTimesFromEventOnset<0.5)
                vectorStrength = np.empty(nRate)
                pValRayleigh = np.empty(nRate)
                for indRate, thisRate in enumerate(possibleRate):
                    # FINISH THIS
                    # - For each rate find spiketimes (use indexLimitsEachTrial)
                    trialsThisRate = trialsEachCond[:,indRate]
                    indexLimitsThisRate = indexLimitsEachTrial[:, trialsThisRate]
                    spikeInds = np.concatenate([np.arange(*pair) for pair in indexLimitsThisRate.T])
                    spiketimesThisRate = spikeTimesFromEventOnset[spikeInds]
                    spiketimesThisRate = spiketimesThisRate[(spiketimesThisRate>period[0]) &
                                                            (spiketimesThisRate<period[1])]
                    if len(spiketimesThisRate):
                        # NOTE: vector strength will be 1 if there is only one spike!
                        strength, phase = signal.vectorstrength(spiketimesThisRate, 1/thisRate)
                        vectorStrength[indRate] = strength
                        radsPerSec = thisRate * 2 * np.pi
                        spikeTimesInRads = (spiketimesThisRate * radsPerSec) % (2 * np.pi)
                        rVal, pVal = studyutils.rayleigh_test(spikeTimesInRads)
                        pValRayleigh[indRate] = pVal
                    else:
                        vectorStrength[indRate] = np.nan
                        pValRayleigh[indRate] = np.nan

        correctedAlphaRayleigh = 0.05/len(possibleRate)
        syncRateInds = np.flatnonzero(pValRayleigh<correctedAlphaRayleigh)
        if len(syncRateInds):
            maxSyncRate = possibleRate[syncRateInds[-1]]
        else:
            maxSyncRate = 0
        celldb.at[indRow, 'amSelectivityPvalOnset'] = amSelectivityPvalOnset[0]
        celldb.at[indRow, 'amSelectivityPvalSustain'] = amSelectivityPvalOnset[1]
        celldb.at[indRow, 'amMaxSyncRate'] = maxSyncRate

        if indCell % 20 == 0:
            print(f'{indCell}/{len(celldbResp)}')
        #print(f'[{indRow}] {str(oneCell)}')

        if 0:
            vsValsStr = ' '.join([f'{vs:0.4f}' for vs in vectorStrength])
            print(vsValsStr)
            pValsStr = ' '.join([f'{p:0.4f}' for p in pValRayleigh])
            print(pValsStr)
            print()
            plt.clf()
            fRaster = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                             timeRange, trialsEachCond)
            fontWeight = 'bold' #if responsiveThisCell else None
            thisTitle = plt.title(f'[{indRow}]', fontweight=fontWeight)
            plt.show()
            #plt.waitforbuttonpress()
            sys.exit()
            #plt.pause(0.1)

    celldatabase.save_hdf(celldb, newdbPath)
