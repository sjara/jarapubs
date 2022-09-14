"""
Calculate laser response to see if cells are D1 or nD1.
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


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning_latency.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

SAVE_DATA = 1
newdbPath = '/tmp/astrpi_am_tuning.h5'
extraDataFile = '/tmp/am_tuning.npz'
scriptFullPath = os.path.realpath(__file__)

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
#celldbResp = celldbResp.loc[1949:] # Fig
#celldbResp = celldbResp.loc[3094:] # Fig
#celldbResp = celldbResp.loc[3420:] # Fig 
#celldbResp = celldbResp.loc[3124:] # Fig 

celldb['amSelectivityPvalOnset'] = np.nan
celldb['amSelectivityPvalSustain'] = np.nan
celldb['amMaxSyncRate'] = np.nan
cellIndexes = celldbResp.index
amFiringRateEachRateSustain = np.full((len(celldbResp), N_RATE), np.nan)
amVectorStrengthEachRateSustain = np.full((len(celldbResp), N_RATE), np.nan)

indCell = -1
for indRow, dbRow in celldbResp.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('am')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    detectorOnsetTimes = ephysData['events']['soundDetectorOn']
    #eventOnsetTimes = spikesanalysis.minimum_event_onset_diff(detectorOnsetTimes, 0.5)
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

        # -- For the sustained period, get average rate and synchronization --
        if indPeriod==1:
            avgFiringEachRate = np.empty(nRate)
            vectorStrength = np.empty(nRate)
            pValRayleigh = np.empty(nRate)
            for indRate, thisRate in enumerate(possibleRate):
                # -- Estimate average firing for each rate --
                trialsThisRate = trialsEachCond[:,indRate]
                spikeCountThisRate = spikeCountMat[trialsThisRate]
                avgFiringEachRate[indRate] = np.mean(spikeCountThisRate)/periodDuration[indPeriod]

                # - - For each rate find spiketimes (use indexLimitsEachTrial) --
                indexLimitsThisRate = indexLimitsEachTrial[:, trialsThisRate]
                spikeInds = np.concatenate([np.arange(*pair) for pair in indexLimitsThisRate.T])
                spiketimesThisRate = spikeTimesFromEventOnset[spikeInds]
                spiketimesThisRate = spiketimesThisRate[(spiketimesThisRate>period[0]) &
                                                        (spiketimesThisRate<period[1])]
                # -- Evaluate synchronization (during sustained) --
                if len(spiketimesThisRate)>1:
                    # NOTE: vector strength is set to NaN if there is only one spike.
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
    amFiringRateEachRateSustain[indCell, :] = avgFiringEachRate
    amVectorStrengthEachRateSustain[indCell, :] = vectorStrength

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
        import matplotlib.gridspec as gridspec
        gsMain = gridspec.GridSpec(1, 3, width_ratios = [0.4, 0.2, 0.2])
        plt.subplot(gsMain[0]) 
        (pRaster,hcond,zline) = extraplots.raster_plot(spikeTimesFromEventOnset,
                                                       indexLimitsEachTrial,
                                                       timeRange, trialsEachCond)
        plt.setp(pRaster, ms=2)
        fontWeight = 'bold' #if responsiveThisCell else None
        thisTitle = plt.title(f'[{indRow}]', fontweight=fontWeight)
        yVals = np.log2(possibleRate)
        plt.subplot(gsMain[1])
        plt.plot(avgFiringEachRate, np.log2(possibleRate), '-o', color='k')
        plt.xlabel('Firing rate (spk/s)')
        plt.ylim([yVals[0]-0.2, yVals[-1]+0.2])
        plt.subplot(gsMain[2])
        plt.plot(vectorStrength, np.log2(possibleRate), '-o', color='k')
        plt.xlabel('Vector strength')
        plt.xlim([0,1])
        plt.ylim([yVals[0]-0.2, yVals[-1]+0.2])
        plt.show()
        #plt.waitforbuttonpress()
        sys.exit()
        #plt.pause(0.1)

if SAVE_DATA:
    celldatabase.save_hdf(celldb, newdbPath)
    np.savez(extraDataFile, script=scriptFullPath, cellIndexes= cellIndexes,
             possibleRate = possibleRate,
             amFiringRateEachRateSustain=amFiringRateEachRateSustain,
             amVectorStrengthEachRateSustain=amVectorStrengthEachRateSustain)
    print(f'Saved {extraDataFile}')


'''
nToneResponsiveD1 = np.sum((celldb.toneMinPval<(0.05/16)) & (celldb.laserPvalB200R50<(0.01)))   
nToneResponsiveND1 = np.sum((celldb.toneMinPval<(0.05/16)) & (celldb.laserPvalB200R50>(0.05)))   

nAMResponsiveD1 = np.sum((celldb.amMinPval<(0.05/11)) & (celldb.laserPvalB200R50<(0.01)))   
nAMResponsiveND1 = np.sum((celldb.amMinPval<(0.05/11)) & (celldb.laserPvalB200R50>(0.05)))   

nToneSelectiveD1 = np.sum((celldb.toneSelectivityPval<0.01) & (celldb.laserPvalB200R50<(0.01)))
nToneSelectiveND1 = np.sum((celldb.toneSelectivityPval<0.01) & (celldb.laserPvalB200R50>(0.05)))

nAMSelectiveD1 = np.sum((celldb.amSelectivityPval<0.01) & (celldb.laserPvalB200R50<(0.01)))
nAMSelectiveND1 = np.sum((celldb.amSelectivityPval<0.01) & (celldb.laserPvalB200R50>(0.05)))  

100*nToneSelectiveD1/nToneResponsiveD1
100*nToneSelectiveND1/nToneResponsiveND1

100*nAMSelectiveD1/nAMResponsiveD1
100*nAMSelectiveND1/nAMResponsiveND1


'''
