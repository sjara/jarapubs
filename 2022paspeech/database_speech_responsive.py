"""
Calculate responsiveness to speech sounds (FT/VOT).
"""
import os
import sys
import studyparams
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

from importlib import reload
reload(extraplots)

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)

allSubjects = studyparams.EPHYS_MICE
#subject = 'feat004'


for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    dbPath = os.path.join(databaseDir, f'{subject}_paspeech_am_tuning.h5')
    celldb = celldatabase.load_hdf(dbPath)
    nCells = len(celldb)

    newdbPath = os.path.join('/tmp', f'{subject}_paspeech_speech_pval.h5')

    periodsName = ['base200', 'respOnset', 'respSustained']
    allPeriods = [ [-0.2, 0], [0, 0.12] , [0.12, 0.24] ]
    periodDuration = [x[1]-x[0] for x in allPeriods]
    N_FT = 4 # HARDCODED
    N_VOT = 4 #HARDCODED

    pValsEachCellOnsetFT = np.empty((nCells, N_FT))
    pValsEachCellSustainFT = np.empty((nCells, N_FT))
    minPvalEachCellOnsetFT = np.full(nCells, np.nan)
    minPvalEachCellSustainFT = np.full(nCells, np.nan)
    minPvalIndexEachCellOnsetFT = np.full(nCells, -1)
    minPvalIndexEachCellSustainFT = np.full(nCells, -1)

    pValsEachCellOnsetVOT = np.empty((nCells, N_VOT))
    pValsEachCellSustainVOT = np.empty((nCells, N_VOT))
    minPvalEachCellOnsetVOT = np.full(nCells, np.nan)
    minPvalEachCellSustainVOT = np.full(nCells, np.nan)
    minPvalIndexEachCellOnsetVOT = np.full(nCells, -1)
    minPvalIndexEachCellSustainVOT = np.full(nCells, -1)

    firingRateEachCellBase = np.full(nCells, np.nan)
    bestFiringRateEachCellOnsetFT = np.full(nCells, np.nan)
    bestFiringRateEachCellSustainFT = np.full(nCells, np.nan)
    bestIndexEachCellOnsetFT = np.full(nCells, -1)
    bestIndexEachCellSustainFT = np.full(nCells, -1)
    maxFiringRateEachCellOnsetFT = np.full(nCells, np.nan)
    minFiringRateEachCellOnsetFT = np.full(nCells, np.nan)
    maxFiringRateEachCellSustainFT = np.full(nCells, np.nan)
    minFiringRateEachCellSustainFT = np.full(nCells, np.nan)

    bestFiringRateEachCellOnsetVOT = np.full(nCells, np.nan)
    bestFiringRateEachCellSustainVOT = np.full(nCells, np.nan)
    bestIndexEachCellOnsetVOT = np.full(nCells, -1)
    bestIndexEachCellSustainVOT = np.full(nCells, -1)
    maxFiringRateEachCellOnsetVOT = np.full(nCells, np.nan)
    minFiringRateEachCellOnsetVOT = np.full(nCells, np.nan)
    maxFiringRateEachCellSustainVOT = np.full(nCells, np.nan)
    minFiringRateEachCellSustainVOT = np.full(nCells, np.nan)

    indCell = -1
    for indRow, dbRow in celldb.iterrows():
        indCell += 1
        oneCell = ephyscore.Cell(dbRow)

        if 'FTVOTBorders' not in oneCell.dbRow.sessionType:
            print(f'[{indRow}] does not have a FT-VOT session')
            continue

        ephysData, bdata = oneCell.load('FTVOTBorders')

        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        timeRange = [-0.4, 0.55]  # In seconds

        # -- Test if trials from behavior don't match ephys --
        ''' shouldn't need this for speech
        if (len(rateEachTrial) > len(eventOnsetTimes)) or \
           (len(rateEachTrial) < len(eventOnsetTimes)-1):
            print(f'[{indRow}] Warning! BevahTrials ({len(rateEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')
            continue
        if len(rateEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(rateEachTrial)]
        '''

        FTParamsEachTrial = bdata['targetFTpercent']
        possibleFTParams = np.unique(FTParamsEachTrial)
        nFT = len(possibleFTParams)
        trialsEachFTCond = behavioranalysis.find_trials_each_type(FTParamsEachTrial, possibleFTParams)

        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        nVOT = len(possibleVOTParams)
        trialsEachVOTCond = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)


        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        meanFiringEachPeriod = np.empty(len(allPeriods))
        spikesEachTrialEachPeriod = []
        for indPeriod, period in enumerate(allPeriods):
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            spikesEachTrialEachPeriod.append(spikesEachTrial)

        firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

        # Calculate mean firing rates and responsiveness for each FT
        meanFiringRateBaseFT = np.empty(nFT)
        meanFiringRateOnsetFT = np.empty(nFT)
        pValEachCondOnsetFT = np.empty(nFT)
        meanFiringRateSustainFT = np.empty(nFT)
        pValEachCondSustainFT = np.empty(nFT)
        for indcond, thisCond in enumerate(possibleFTParams):
            trialsThisCond = trialsEachFTCond[:,indcond]
            firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
            firingRateOnset = spikesEachTrialEachPeriod[1][trialsThisCond]/periodDuration[1]
            firingRateSustain = spikesEachTrialEachPeriod[2][trialsThisCond]/periodDuration[2]
            try:
                wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateOnset)
            except ValueError:
                pValThisCond = 1
            pValEachCondOnsetFT[indcond] = pValThisCond
            try:
                wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateSustain)
            except ValueError:
                pValThisCond = 1
            pValEachCondSustainFT[indcond] = pValThisCond
            meanFiringRateOnsetFT[indcond] = firingRateOnset.mean()
            meanFiringRateSustainFT[indcond] = firingRateSustain.mean()

        indMinPvalOnset = np.argmin(pValEachCondOnsetFT)
        minPvalEachCellOnsetFT[indCell] = pValEachCondOnsetFT[indMinPvalOnset]
        minPvalIndexEachCellOnsetFT[indCell] = indMinPvalOnset
        indBestOnset = np.argmax(np.abs(meanFiringRateOnsetFT-firingRateEachCellBase[indCell]))
        bestFiringRateEachCellOnsetFT[indCell] = meanFiringRateOnsetFT[indBestOnset]
        bestIndexEachCellOnsetFT[indCell] = indBestOnset
        maxFiringRateEachCellOnsetFT[indCell] = np.max(meanFiringRateOnsetFT)
        minFiringRateEachCellOnsetFT[indCell] = np.min(meanFiringRateOnsetFT)

        indMinPvalSustain = np.argmin(pValEachCondSustainFT)
        minPvalEachCellSustainFT[indCell] = pValEachCondSustainFT[indMinPvalSustain]
        minPvalIndexEachCellSustainFT[indCell] = indMinPvalSustain
        indBestSustain = np.argmax(np.abs(meanFiringRateSustainFT-firingRateEachCellBase[indCell]))
        bestFiringRateEachCellSustainFT[indCell] = meanFiringRateSustainFT[indBestSustain]
        bestIndexEachCellSustainFT[indCell] = indBestSustain
        maxFiringRateEachCellSustainFT[indCell] = np.max(meanFiringRateSustainFT)
        minFiringRateEachCellSustainFT[indCell] = np.min(meanFiringRateSustainFT)


        # Calculate mean firing rates and responsiveness for each VOT
        meanFiringRateBaseVOT = np.empty(nVOT)
        meanFiringRateOnsetVOT = np.empty(nVOT)
        pValEachCondOnsetVOT = np.empty(nVOT)
        meanFiringRateSustainVOT = np.empty(nVOT)
        pValEachCondSustainVOT = np.empty(nVOT)
        for indcond, thisCond in enumerate(possibleVOTParams):
            trialsThisCond = trialsEachVOTCond[:,indcond]
            firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
            firingRateOnset = spikesEachTrialEachPeriod[1][trialsThisCond]/periodDuration[1]
            firingRateSustain = spikesEachTrialEachPeriod[2][trialsThisCond]/periodDuration[2]
            try:
                wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateOnset)
            except ValueError:
                pValThisCond = 1
            pValEachCondOnsetVOT[indcond] = pValThisCond
            try:
                wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateSustain)
            except ValueError:
                pValThisCond = 1
            pValEachCondSustainVOT[indcond] = pValThisCond
            meanFiringRateOnsetVOT[indcond] = firingRateOnset.mean()
            meanFiringRateSustainVOT[indcond] = firingRateSustain.mean()

        indMinPvalOnset = np.argmin(pValEachCondOnsetVOT)
        minPvalEachCellOnsetVOT[indCell] = pValEachCondOnsetVOT[indMinPvalOnset]
        minPvalIndexEachCellOnsetVOT[indCell] = indMinPvalOnset
        indBestOnset = np.argmax(np.abs(meanFiringRateOnsetVOT-firingRateEachCellBase[indCell]))
        bestFiringRateEachCellOnsetVOT[indCell] = meanFiringRateOnsetVOT[indBestOnset]
        bestIndexEachCellOnsetVOT[indCell] = indBestOnset
        maxFiringRateEachCellOnsetVOT[indCell] = np.max(meanFiringRateOnsetVOT)
        minFiringRateEachCellOnsetVOT[indCell] = np.min(meanFiringRateOnsetVOT)

        indMinPvalSustain = np.argmin(pValEachCondSustainVOT)
        minPvalEachCellSustainVOT[indCell] = pValEachCondSustainVOT[indMinPvalSustain]
        minPvalIndexEachCellSustainVOT[indCell] = indMinPvalSustain
        indBestSustain = np.argmax(np.abs(meanFiringRateSustainVOT-firingRateEachCellBase[indCell]))
        bestFiringRateEachCellSustainVOT[indCell] = meanFiringRateSustainVOT[indBestSustain]
        bestIndexEachCellSustainVOT[indCell] = indBestSustain
        maxFiringRateEachCellSustainVOT[indCell] = np.max(meanFiringRateSustainVOT)
        minFiringRateEachCellSustainVOT[indCell] = np.min(meanFiringRateSustainVOT)




    celldb['FTMinPvalOnset'] = minPvalEachCellOnsetFT
    celldb['FTIndexMinPvalOnset'] = minPvalIndexEachCellOnsetFT
    celldb['FTMinPvalSustain'] = minPvalEachCellSustainFT
    celldb['FTIndexMinPvalOnset'] = minPvalIndexEachCellSustainFT
    celldb['VOTMinPvalOnset'] = minPvalEachCellOnsetVOT
    celldb['VOTIndexMinPvalOnset'] = minPvalIndexEachCellOnsetVOT
    celldb['VOTMinPvalSustain'] = minPvalEachCellSustainVOT
    celldb['VOTIndexMinPvalOnset'] = minPvalIndexEachCellSustainVOT

    celldb['speechFiringRateBaseline'] = firingRateEachCellBase
    celldb['FTFiringRateBestOnset'] = bestFiringRateEachCellOnsetFT
    celldb['FTIndexBestOnset'] = bestIndexEachCellOnsetFT
    celldb['FTFiringRateBestSustain'] = bestFiringRateEachCellSustainFT
    celldb['FTIndexBestSustain'] = bestIndexEachCellSustainFT
    celldb['VOTFiringRateBestOnset'] = bestFiringRateEachCellOnsetVOT
    celldb['VOTIndexBestOnset'] = bestIndexEachCellOnsetVOT
    celldb['VOTFiringRateBestSustain'] = bestFiringRateEachCellSustainVOT
    celldb['VOTIndexBestSustain'] = bestIndexEachCellSustainVOT

    celldb['FTFiringRateMaxOnset'] = maxFiringRateEachCellOnsetFT
    celldb['FTFiringRateMinOnset'] = minFiringRateEachCellOnsetFT
    celldb['FtFiringRateMaxSustain'] = maxFiringRateEachCellSustainFT
    celldb['FTFiringRateMinSustain'] = minFiringRateEachCellSustainFT
    celldb['VOTFiringRateMaxOnset'] = maxFiringRateEachCellOnsetVOT
    celldb['VOTFiringRateMinOnset'] = minFiringRateEachCellOnsetVOT
    celldb['VOTFiringRateMaxSustain'] = maxFiringRateEachCellSustainVOT
    celldb['VOTFiringRateMinSustain'] = minFiringRateEachCellSustainVOT

    celldatabase.save_hdf(celldb, newdbPath)
