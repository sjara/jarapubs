"""
Calculate frequency tuning measurements.
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
from scipy import optimize
from scipy import signal

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

import studyutils

from importlib import reload
reload(studyutils)
reload(spikesanalysis)


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

newdbPath = '/tmp/astrpi_freq_tuning_latency.h5'

N_FREQ = 16 # HARDCODED
correctedAlpha = 0.05/N_FREQ
celldbResp = celldb[celldb.toneMinPval < correctedAlpha]

#celldbResp = celldbResp.loc[1391:]  # Example in figure
#celldbResp = celldbResp.loc[1999:]  # Potential Example in figure
#celldbResp = celldbResp.loc[3068:]  # Good response but latency problems
#celldbResp = celldbResp.loc[739:]  # Good response but very short latency
#celldbResp = celldbResp.loc[1381:]  # Supressed response
#celldbResp = celldbResp.loc[1728:]  # Supressed response

indCell = -1
for indRow, dbRow in celldbResp.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('tuningCurve')

    spikeTimes = ephysData['spikeTimes']
    #eventOnsetTimes = ephysData['events']['stimOn']
    detectorOnsetTimes = ephysData['events']['soundDetectorOn']
    eventOnsetTimes = spikesanalysis.minimum_event_onset_diff(detectorOnsetTimes, 0.5)
    latencyTimeRange = [-0.2, 0.2]

    freqEachTrial = bdata['currentFreq']
    intensityEachTrial = bdata['currentIntensity']
    # Remove last stim from ephys if not saved in behavior file
    if len(freqEachTrial) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(freqEachTrial)]

    if 0:
        # -- OPTION 1: best freq, above X dB --
        bestFreqInd = dbRow['toneIndexBest']
        thisFreq = np.unique(bdata['currentFreq'])[bestFreqInd]
        useIntensityAbove = 40
        selectedTrials = (bdata['currentFreq']==thisFreq) & (bdata['currentIntensity']>useIntensityAbove)
    else:
        # -- OPTION 2: freq+intensity within FRA
        # -- Calculate average firing for each freq-intensity combination
        responsePeriod = [0, 0.1]
        respPeriodDuration = responsePeriod[1]-responsePeriod[0]
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, latencyTimeRange)
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 responsePeriod)
        nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it
        possibleFreq = np.unique(freqEachTrial)
        nFreq = len(possibleFreq)
        possibleIntensity = np.unique(intensityEachTrial)
        nIntensity = len(possibleIntensity)
        trialsEachComb = behavioranalysis.find_trials_each_combination(freqEachTrial,
                                                                       possibleFreq,
                                                                       intensityEachTrial,
                                                                       possibleIntensity)
        # -- Calculate average firing for each freq-intensity combination
        firingRateRespMap = np.empty((nIntensity, nFreq))
        for indFreq, frequency in enumerate(possibleFreq):
            for indInt, intensity in enumerate(possibleIntensity):
                trialsThisComb = trialsEachComb[:, indFreq, indInt]
                nSpikesThisComb = nSpikesEachTrial[trialsThisComb]
                firingRateRespMap[indInt, indFreq] = np.mean(nSpikesThisComb)/respPeriodDuration
                
        firingRateBaseline = dbRow['toneFiringRateBaseline']
        if np.mean(firingRateRespMap) > firingRateBaseline:
            fra, respThreshold = studyutils.calculate_fra(firingRateRespMap, firingRateBaseline)
            invertedResponse = False
        else:
            invertedResponse = True
            print(f'Inverting FRA for cell {indRow}')
            invertedMap = firingRateRespMap.max()-firingRateRespMap
            invertedBase = firingRateRespMap.max()-firingRateBaseline
            fra, respThreshold = studyutils.calculate_fra(invertedMap, invertedBase)
        selectedTrials = trialsEachComb[:,fra.T].sum(axis=1).astype(bool)

        #meshF, meshI = np.meshgrid(possibleFreq, possibleIntensity)
        #flatF = meshF[fra]
        #flatI = meshI[fra]
        #selectedTrials = np.zeros(len(freqEachTrial), dtype=bool)
        #for oneF, oneI in zip(flatF, flatI):
        #    trialsThisComb
        #    selectedTrials = (selectedTrials|

    indexLimitsEachTrial = indexLimitsEachTrial[:, selectedTrials]
    
    # -- Estimate latency --
    smoothWin = signal.windows.hann(11)
    respLatency, interim = spikesanalysis.response_latency(spikeTimesFromEventOnset,
                                                           indexLimitsEachTrial, latencyTimeRange,
                                                           win=smoothWin, invert=invertedResponse)
    #print(f'[{indRow}]: {1e3*respLatency:0.1f} ms')

    # -- Estimate onset vs sustain response --
    onsetSustainRanges = [0, 0.05, 0.1]
    onsetDuration = onsetSustainRanges[1] - onsetSustainRanges[0]
    sustainDuration = onsetSustainRanges[2] - onsetSustainRanges[1]
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                             indexLimitsEachTrial,
                                                             onsetSustainRanges)
    firingRateOnset = np.mean(spikeCountMat[:,0])/onsetDuration
    firingRateSustain = np.mean(spikeCountMat[:,1])/sustainDuration
    
    
    celldb.at[indRow, 'toneLatencySoundDetector'] = respLatency
    celldb.at[indRow, 'toneFiringRateBestOnset'] = firingRateOnset
    celldb.at[indRow, 'toneFiringRateBestSustain'] = firingRateSustain
    
    if indCell % 10 == 0:
        print(f'{indCell}/{len(celldbResp)}')
        pass
    #print(f'[{indRow}] {str(oneCell)}')
    
    if 0:
        plt.clf()
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials],
                                                  latencyTimeRange)
        markerSize = 2
        ax0 = plt.subplot(2,1,1)
        plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=markerSize)
        plt.axvline(0, color='b')
        plt.axvline(respLatency, color='r')
        plt.xlim([-0.2, 0.2])
        plt.title(f'[{indRow}]: {1e3*respLatency:0.1f} ms')
        ax1 = plt.subplot(2,1,2, sharex=ax0)
        plt.plot(interim['timeVec'], interim['psth'])
        plt.show()
        #sys.exit()
        plt.waitforbuttonpress()
        #sys.exit()
        #plt.pause(0.5)

celldatabase.save_hdf(celldb, newdbPath)
 
'''
clf(); hist(celldb.toneLatencyBehav[~np.isnan(celldb.toneLatencyBehav)],200);
clf(); hist(celldb.toneLatencySoundDetector[~np.isnan(celldb.toneLatencySoundDetector)],200); 
'''
