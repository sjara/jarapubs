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

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

import studyutils

from importlib import reload
reload(studyutils)


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_am_pval.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

newdbPath = '/tmp/astrpi_freq_tuning.h5'

#periodsName = ['base200', 'resp100']
#allPeriods = [ [-0.2, 0], [0, 0.1] ]
#periodDuration = [x[1]-x[0] for x in allPeriods]
#meanFiringEachPeriodEachCell = np.empty((nCells, len(allPeriods)))

responsePeriod = [0, 0.1]
respPeriodDuration = responsePeriod[1]-responsePeriod[0]

N_FREQ = 16 # HARDCODED
correctedAlpha = 0.05/N_FREQ

'''
toneSelectivityPvalEachCell = np.full(nCells, np.nan)
toneIntensityThresholdEachCell = np.full(nCells, np.nan)
toneCharactFreqEachCell = np.full(nCells, np.nan)
toneGaussianRsquare = np.full(nCells, np.nan)
toneGaussianSigma = np.full(nCells, np.nan)
toneGaussianX0 = np.full(nCells, np.nan)
toneGaussianA = np.full(nCells, np.nan)
toneGaussianY0 = np.full(nCells, np.nan)
'''

celldbResp = celldb[celldb.toneMinPval < correctedAlpha]

celldb['toneSelectivityPval'] = np.nan
celldb['toneFiringRateBestOnset'] = np.nan
celldb['toneFiringRateBestSustain'] = np.nan
celldb['toneCharactFreqEachCell'] = np.nan 
celldb['toneGaussianRsquare'] = np.nan
celldb['toneGaussianSigma'] = np.nan
celldb['toneGaussianX0'] = np.nan
celldb['toneGaussianA'] = np.nan
celldb['toneGaussianY0'] = np.nan

#celldbResp = celldbResp.loc[2858:]
#celldbResp = celldbResp.loc[516:]
#celldbResp = celldbResp.loc[650:]
#celldbResp = celldbResp.loc[739:]
#celldbResp = celldbResp.loc[926:]  # Gives an error
#celldbResp = celldbResp.loc[1317:]  # Good test
#celldbResp = celldbResp.loc[1381:]  # Inverted test
#celldbResp = celldbResp.loc[717:]  # Gives an error
#celldbResp = celldbResp.loc[710:]  # Example in figure
#celldbResp = celldbResp.loc[1317:]  # Potential Example in figure
#celldbResp = celldbResp.loc[1391:]  # Example in figure
#celldbResp = celldbResp.loc[1999:]  # Potential Example in figure

indCell = -1
for indRow, dbRow in celldbResp.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('tuningCurve')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.3, 0.6]  # In seconds

    freqEachTrial = bdata['currentFreq']
    intensityEachTrial = bdata['currentIntensity']
    # Remove last stim from ephys if not saved in behavior file
    if len(freqEachTrial) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(freqEachTrial)]
    
    possibleFreq = np.unique(freqEachTrial)
    nFreq = len(possibleFreq)
    possibleIntensity = np.unique(intensityEachTrial)
    nIntensity = len(possibleIntensity)
    
    trialsEachComb = behavioranalysis.find_trials_each_combination(freqEachTrial,
                                                                   possibleFreq,
                                                                   intensityEachTrial,
                                                                   possibleIntensity)

    #trialsEachCond = behavioranalysis.find_trials_each_type(freqEachTrial, possibleFreq)
    trialsEachCond = trialsEachComb.sum(axis=2).astype('bool')

    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                             indexLimitsEachTrial,
                                                             responsePeriod)
    nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

    # -- Calculate nSpikes for each freq to test selectivity --
    nSpikesEachFreq = []
    for indFreq, frequency in enumerate(possibleFreq):
        nSpikesEachFreq.append(nSpikesEachTrial[trialsEachCond[:,indFreq]])
    kStat, pValKruskal = stats.kruskal(*nSpikesEachFreq)
    #print(pValKruskal)
    #celldb.at[indRow, 'toneSelectivityPval'] = pValKruskal
    #toneSelectivityPvalEachCell[indRow] = pValKruskal
 
    # -- Calculate average firing for each freq-intensity combination
    firingRateRespMap = np.empty((nIntensity, nFreq))
    for indFreq, frequency in enumerate(possibleFreq):
        for indInt, intensity in enumerate(possibleIntensity):
            trialsThisComb = trialsEachComb[:, indFreq, indInt]
            nSpikesThisComb = nSpikesEachTrial[trialsThisComb]
            firingRateRespMap[indInt, indFreq] = np.mean(nSpikesThisComb)/respPeriodDuration
            #nTrialsThisComb = len(nSpikesThisComb)

    firingRateBaseline = dbRow['toneFiringRateBaseline']
    if np.mean(firingRateRespMap) > firingRateBaseline:
        fra, respThreshold = studyutils.calculate_fra(firingRateRespMap, firingRateBaseline)
    else:
        print(f'Inverting FRA for cell {indRow}')
        invertedMap = firingRateRespMap.max()-firingRateRespMap
        invertedBase = firingRateRespMap.max()-firingRateBaseline
        fra, respThreshold = studyutils.calculate_fra(invertedMap, invertedBase)
    intenThresholdInd, cfInd = studyutils.calculate_intensity_threshold(fra)
    #print(intenThresholdInd, cfInd)
    #toneIntensityThresholdEachCell[indRow] = intenThresholdInd
    #toneCharactFreqEachCell[indRow] = cfInd
    if intenThresholdInd is not None:
        maxIntensityInd = np.argmax(firingRateRespMap.sum(axis=1))
        (lowSlope, highSlope) = studyutils.calculate_fra_slopes(fra, possibleIntensity,
                                                                possibleFreq, cfInd,
                                                                intenThresholdInd, maxIntensityInd)

    '''
    # -- Fit a Gaussian to a specific intensity --
    #selTntensity = possibleIntensity[intenThresholdInd]+10
    selTntensity = 60
    indSelIntensity = (np.abs(possibleIntensity - selTntensity)).argmin()
    nTrialsSelIntensity = trialsEachComb[:, :, indSelIntensity].sum() 
    nSpikesSelIntensity = np.full(nTrialsSelIntensity, np.nan)
    logFreqsSelIntensity = np.full(nTrialsSelIntensity, np.nan)
    possibleLogFreq = np.log2(possibleFreq)
    trialIndex = 0
    for indFreq, frequency in enumerate(possibleFreq):
        trialsThisComb = trialsEachComb[:, indFreq, indSelIntensity]
        nTrialsThisComb = trialsThisComb.sum()
        trialRange = np.arange(nTrialsThisComb) + trialIndex
        nSpikesSelIntensity[trialRange] = nSpikesEachTrial[trialsThisComb]
        logFreqsSelIntensity[trialRange] = possibleLogFreq[indFreq]
        trialIndex += nTrialsThisComb
    '''
    logFreqsSelIntensity = np.log2(freqEachTrial)
    possibleLogFreq = np.log2(possibleFreq)
    nSpikesSelIntensity = nSpikesEachTrial
    # PARAMS: a, x0, sigma, y0
    if dbRow.toneFiringRateBest >= dbRow.toneFiringRateBaseline:
        p0 = [1, possibleLogFreq[N_FREQ//2], 1, firingRateBaseline]
        bounds = ([0, possibleLogFreq[0], 0, 0],
                  [np.inf, possibleLogFreq[-1], np.inf, np.inf])
    else:
        p0 = [-1, possibleLogFreq[N_FREQ//2], 1, firingRateBaseline]
        bounds = ([-np.inf, possibleLogFreq[0], 0, 0],
                  [0, possibleLogFreq[-1], np.inf, np.inf])
    try:
        popt, pcov = optimize.curve_fit(studyutils.gaussian, logFreqsSelIntensity,
                                        nSpikesSelIntensity, p0=p0, bounds=bounds)
    except RuntimeError:
        print("Could not fit gaussian curve to tuning data.")
        popt = None
        Rsquared = None
    else:
        gaussianResp = studyutils.gaussian(logFreqsSelIntensity, *popt)
        residuals = nSpikesSelIntensity - gaussianResp
        ssquared = np.sum(residuals**2)
        ssTotal = np.sum((nSpikesSelIntensity-np.mean(nSpikesSelIntensity))**2)
        Rsquared = 1 - (ssquared/ssTotal)
        fullWidthHalfMax = 2.355*popt[2] # Sigma is popt[2]
    
    #print(Rsquared)
    #print(fullWidthHalfMax) # Sigma
    '''
    toneGaussianRsquare[indRow] = Rsquared
    toneGaussianA[indRow] = popt[0]
    toneGaussianX0[indRow] = popt[1]
    toneGaussianSigma[indRow] = popt[2]
    toneGaussianY0[indRow] = popt[3]
    '''

    # -- Estimate BW10 --
    RsquaredBW10 = np.nan
    fullWHMaxBW10 = np.nan
    RsquaredBW40 = np.nan
    fullWHMaxBW40 = np.nan
    if (intenThresholdInd is not None) and (np.mean(firingRateRespMap) > firingRateBaseline):
        intensityThreshold = possibleIntensity[intenThresholdInd]
        if intensityThreshold<65:
            dBabove = 10
            (poptBW10, RsquaredBW10, fullWHMaxBW10) = studyutils.calculate_BWx(intensityEachTrial,
                                                                               freqEachTrial,
                                                                               nSpikesEachTrial,
                                                                               intensityThreshold,
                                                                               dBabove)
        if intensityThreshold<35:
            dBabove = 40
            (poptBW40, RsquaredBW40, fullWHMaxBW40) = studyutils.calculate_BWx(intensityEachTrial,
                                                                               freqEachTrial,
                                                                               nSpikesEachTrial,
                                                                               intensityThreshold,
                                                                               dBabove)
    # -- Estimate latency (from taskontrol timing) --
    bestFreqInd = dbRow['toneIndexBest']
    thisFreq = np.unique(bdata['currentFreq'])[bestFreqInd]
    useIntensityAbove = 40
    selectedTrials = (bdata['currentFreq']==thisFreq) & (bdata['currentIntensity']>useIntensityAbove)
    latencyTimeRange = [-0.2, 0.2]
    indexLimitsSelectedTrials = indexLimitsEachTrial[:,selectedTrials]
    smoothWin = signal.windows.hann(11)
    respLatency, interim = spikesanalysis.response_latency(spikeTimesFromEventOnset,
                                                           indexLimitsSelectedTrials, latencyTimeRange,
                                                           win=smoothWin)
    #print(f"{respLatency:0.4f}")

    '''
    # -- Estimate onset vs sustain response  (this is not in database_tone_latency.py)--
    onsetSustainRanges = [0, 0.05, 0.1]
    #onsetSustainRanges = [0.01, 0.05, 0.09]
    onsetDuration = onsetSustainRanges[1] - onsetSustainRanges[0]
    sustainDuration = onsetSustainRanges[2] - onsetSustainRanges[1]
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                             indexLimitsEachTrial,
                                                             onsetSustainRanges)
    bestFreqInd = dbRow['toneIndexBest']
    trialsThisFreq = trialsEachCond[:,bestFreqInd]
    firingRateOnset = np.mean(spikeCountMat[trialsThisFreq, 0])/onsetDuration
    firingRateSustain = np.mean(spikeCountMat[trialsThisFreq, 1])/sustainDuration
    '''
    
    celldb.at[indRow, 'toneSelectivityPval'] = pValKruskal
    #celldb.at[indRow, 'toneFiringRateBestOnset'] = firingRateOnset
    #celldb.at[indRow, 'toneFiringRateBestSustain'] = firingRateSustain
    celldb.at[indRow, 'toneBW10'] = fullWHMaxBW10
    celldb.at[indRow, 'toneBW10Rsquared'] = RsquaredBW10
    celldb.at[indRow, 'toneBW40'] = fullWHMaxBW40
    celldb.at[indRow, 'toneBW40Rsquared'] = RsquaredBW40
    celldb.at[indRow, 'toneLatencyBehav'] = respLatency
    if intenThresholdInd is not None:
        celldb.at[indRow, 'toneIntensityThreshold'] = possibleIntensity[intenThresholdInd]
        celldb.at[indRow, 'toneCharactFreq'] = possibleFreq[cfInd]
        celldb.at[indRow, 'toneFRAlowSlope'] = lowSlope
        celldb.at[indRow, 'toneFRAhighSlope'] = highSlope
    if Rsquared is not None:
        celldb.at[indRow, 'toneGaussianRsquare'] = Rsquared
        celldb.at[indRow, 'toneGaussianA'] = popt[0]
        celldb.at[indRow, 'toneGaussianX0'] = popt[1]
        celldb.at[indRow, 'toneGaussianSigma'] = popt[2]
        celldb.at[indRow, 'toneGaussianY0'] = popt[3]
    
    if indCell % 10 == 0:
        print(f'{indCell}/{len(celldbResp)}')
        pass
    #print(f'[{indRow}] {str(oneCell)}')
    
    if 0:
        plt.clf()
        plt.subplot2grid((3,2),[0,0], rowspan=3)
        fRaster = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                         timeRange, trialsEachCond)
        fontWeight = 'bold' #if responsiveThisCell else None
        thisTitle = plt.title(f'[{indRow}]', fontweight=fontWeight)
        axFRA = plt.subplot2grid((3,2),[0,1])
        imExtent = None
        maxFR = None
        plt.imshow(firingRateRespMap, interpolation='nearest',
                   cmap='Blues', origin='lower', extent=imExtent, vmax=maxFR)
        plt.colorbar()
        plt.subplot2grid((3,2),[1,1], sharex=None)
        xvals = np.linspace(possibleLogFreq[0], possibleLogFreq[-1], 60)
        yvals = studyutils.gaussian(xvals, *popt)
        plt.plot(xvals, yvals, '.')
        plt.ylabel(f'{Rsquared:0.4f}')
        plt.subplot2grid((3,2),[2,1], sharex=None)
        plt.imshow(fra.astype(int), interpolation='nearest',
                   cmap='Blues', origin='lower', extent=imExtent)
        plt.colorbar()
        #, extent=[0, nFreq-1, nIntensity-1, 0]
        plt.show()
        #plt.waitforbuttonpress()
        sys.exit()
        #plt.pause(0.1)


'''        
celldb['toneSelectivityPval'] = toneSelectivityPval
celldb['toneIntensityThreshold'] = toneIntensityThresholdEachCell
celldb['toneCharactFreq'] =  toneCharactFreqEachCell
celldb['toneGaussianRsquare'] = toneGaussianRsquare
celldb['toneGaussianSigma'] = toneGaussianSigma
celldb['toneGaussianX0'] = toneGaussianX0
celldb['toneGaussianA'] = toneGaussianA
celldb['toneGaussianY0'] = toneGaussianY0
'''

celldatabase.save_hdf(celldb, newdbPath)


'''
plt.clf()
x = np.arange(-10,10,0.1)
a=2; x0=1; sigma=2; y0=1;
y = a*np.exp(-(x-x0)**2/(2*sigma**2))+y0
plt.plot(x,y,'.');

y = a*np.arange(0.001, 0.999, 0.01)+y0
xsq = -2*(sigma**2) * np.log((y-y0)/a)
x = (x0 + np.sqrt(xsq))
plt.plot(x,y,'.');
'''

'''
goodFit = celldb.toneGaussianRsquare>0.01
bandwidthGoodFit = 2.355*celldb.toneGaussianSigma[goodFit]
plt.clf()
plt.hist(bandwidthGoodFit,60) 

'''
