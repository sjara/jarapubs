"""
Calculate frequency tuning
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

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)

#allSubjects = studyparams.EPHYS_MICE
allSubjects = ['feat009', 'feat010']

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    dbPath = os.path.join(databaseDir, f'{subject}_paspeech_tones_pval.h5')
    celldb = celldatabase.load_hdf(dbPath)
    nCells = len(celldb)

    newdbPath = os.path.join('/tmp', f'{subject}_paspeech_freq_tuning.h5')

    responsePeriod = [0, 0.1]
    respPeriodDuration = responsePeriod[1]-responsePeriod[0]

    N_FREQ = 16 # HARDCODED

    correctedAlpha = 0.05/N_FREQ

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


    indCell = -1
    for indRow, dbRow in celldbResp.iterrows():
        indRow
        #dbRow = celldb.loc[570]
        indCell += 1
        oneCell = ephyscore.Cell(dbRow)

        ephysData, bdata = oneCell.load('pureTones')

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
        fra, respThreshold = studyutils.calculate_fra(firingRateRespMap, firingRateBaseline)
        intenThresholdInd, cfInd = studyutils.calculate_intensity_threshold(fra)
        #print(intenThresholdInd, cfInd)
        #toneIntensityThresholdEachCell[indRow] = intenThresholdInd
        #toneCharactFreqEachCell[indRow] = cfInd

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

        # -- Estimate onset vs sustain response --
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

        celldb.at[indRow, 'toneSelectivityPval'] = pValKruskal
        celldb.at[indRow, 'toneFiringRateBestOnset'] = firingRateOnset
        celldb.at[indRow, 'toneFiringRateBestSustain'] = firingRateSustain
        if intenThresholdInd is not None:
            celldb.at[indRow, 'toneIntensityThreshold'] = possibleIntensity[intenThresholdInd]
            celldb.at[indRow, 'toneCharactFreq'] = possibleFreq[cfInd]
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
            plt.imshow(firingRateRespMap, interpolation='nearest',
                       cmap='Blues', origin='lower', extent=imExtent)
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


    celldatabase.save_hdf(celldb, newdbPath)
