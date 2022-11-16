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

#allSubjects = studyparams.EPHYS_MICE
#allSubjects = ['feat009', 'feat010']
allSubjects = studyparams.TEST_MOUSE


for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    dbPath = os.path.join(databaseDir, f'{subject}_paspeech_am_tuning.h5')
    celldb = celldatabase.load_hdf(dbPath)
    nCells = len(celldb)

    newdbPath = os.path.join(databaseDir, f'{subject}_paspeech_speech_pval.h5')

    periodsName = ['base200', 'respOnset', 'respSustained']
    allPeriods = [ [-0.2, 0], [0, 0.12] , [0.12, 0.24] ] #try with shorter period for onset response. I think this is generous.onset vs. sustain?


    periodDuration = [x[1]-x[0] for x in allPeriods]

    #N_stim = 12
    N_FT = 4 # HARDCODED
    N_VOT = 4 #HARDCODE

    pValsEachCellOnset = np.ones((nCells, N_FT, N_VOT))
    pValsEachCellSustain = np.ones((nCells, N_FT, N_VOT))
    minPvalEachCellOnset = np.full(nCells, np.nan)
    minPvalEachCellSustain = np.full(nCells, np.nan)
    minPvalIndexEachCellOnset = np.full((nCells,2), -1)
    minPvalIndexEachCellSustain = np.full((nCells,2), -1)
    firingRateEachCellBase = np.full(nCells, np.nan)
    bestFiringRateEachCellOnset = np.full(nCells, np.nan)
    bestFiringRateEachCellSustain = np.full(nCells, np.nan)
    bestIndexEachCellOnset = np.full((nCells,2), -1)
    bestIndexEachCellSustain = np.full((nCells,2), -1)
    minFiringRate_FT_VOTmin = np.full(nCells, np.nan)
    maxFiringRate_FT_VOTmin = np.full(nCells, np.nan)
    minFiringRate_FT_VOTmax = np.full(nCells, np.nan)
    maxFiringRate_FT_VOTmax = np.full(nCells, np.nan)
    minFiringRate_VOT_FTmin = np.full(nCells, np.nan)
    maxFiringRate_VOT_FTmin = np.full(nCells, np.nan)
    minFiringRate_VOT_FTmax = np.full(nCells, np.nan)
    maxFiringRate_VOT_FTmax = np.full(nCells, np.nan)
    maxFiringRateEachCellOnset = np.full(nCells, np.nan)
    minFiringRateEachCellOnset = np.full(nCells, np.nan)
    maxFiringRateEachCellSustain = np.full(nCells, np.nan)
    minFiringRateEachCellSustain = np.full(nCells, np.nan)


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

        FTParamsEachTrial = bdata['targetFTpercent']

        # -- Test if trials from behavior don't match ephys --
        if (len(FTParamsEachTrial) > len(eventOnsetTimes)) or \
           (len(FTParamsEachTrial) < len(eventOnsetTimes)-1):
            print(f'[{indRow}] Warning! BevahTrials ({len(rateEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')
            continue
        if len(FTParamsEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(FTParamsEachTrial)]


        possibleFTParams = np.unique(FTParamsEachTrial)
        nFT = len(possibleFTParams)
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        nVOT = len(possibleVOTParams)
        nStim = 12

        trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)


        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        meanFiringEachPeriod = np.empty(len(allPeriods))
        spikesEachTrialEachPeriod = []
        for indPeriod, period in enumerate(allPeriods):
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            spikesEachTrialEachPeriod.append(spikesEachTrial)

        firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

        # Calculate mean firing rates and responsiveness for each speech sound (FT-VOT combination)
        meanFiringRateBase = np.empty([nFT, nVOT])
        meanFiringRateOnset = np.empty([nFT, nVOT])
        pValEachCondOnset = np.ones([nFT, nVOT])
        meanFiringRateSustain = np.empty([nFT, nVOT])
        pValEachCondSustain = np.ones([nFT, nVOT])

        for indFT, thisFT in enumerate(possibleFTParams):
            for indVOT, thisVOT in enumerate(possibleVOTParams):
                trialsThisCond = trialsEachCond[:, indFT, indVOT]
                firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
                firingRateOnset = spikesEachTrialEachPeriod[1][trialsThisCond]/periodDuration[1]
                firingRateSustain = spikesEachTrialEachPeriod[2][trialsThisCond]/periodDuration[2]

                try:
                    wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateOnset)
                except ValueError:
                    pValThisCond = 1
                pValEachCondOnset[indFT, indVOT] = pValThisCond
                try:
                    wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateSustain)
                except ValueError:
                    pValThisCond = 1

                pValEachCondSustain[indFT, indVOT] = pValThisCond
                meanFiringRateOnset[indFT, indVOT] = firingRateOnset.mean()
                meanFiringRateSustain[indFT, indVOT] = firingRateSustain.mean()
                meanFiringRateBase[indFT, indVOT] = firingRateBase.mean()


        indMinPvalOnset = np.unravel_index(np.argmin(pValEachCondOnset), pValEachCondOnset.shape)
        minPvalEachCellOnset[indCell] = pValEachCondOnset[indMinPvalOnset]
        minPvalIndexEachCellOnset[indCell,:] = indMinPvalOnset
        indBestOnset = np.unravel_index(np.nanargmax(np.abs(meanFiringRateOnset-firingRateEachCellBase[indCell])), meanFiringRateOnset.shape)
        bestFiringRateEachCellOnset[indCell] = meanFiringRateOnset[indBestOnset]
        bestIndexEachCellOnset[indCell,:] = indBestOnset
        #maxFiringRateEachCellOnset[indCell] = np.nanmax(meanFiringRateOnset)
        #minFiringRateEachCellOnset[indCell] = np.nanmin(meanFiringRateOnset)
        minFiringRate_FT_VOTmin[indCell] = np.nanmin(meanFiringRateOnset[:,0])
        maxFiringRate_FT_VOTmin[indCell] = np.nanmax(meanFiringRateOnset[:,0])
        minFiringRate_FT_VOTmax[indCell] = np.nanmin(meanFiringRateOnset[:,3])
        maxFiringRate_FT_VOTmax[indCell] = np.nanmax(meanFiringRateOnset[:,3])
        minFiringRate_VOT_FTmin[indCell] = np.nanmin(meanFiringRateOnset[0,:])
        maxFiringRate_VOT_FTmin[indCell] = np.nanmax(meanFiringRateOnset[0,:])
        minFiringRate_VOT_FTmax[indCell] = np.nanmin(meanFiringRateOnset[3,:])
        maxFiringRate_VOT_FTmax[indCell] = np.nanmax(meanFiringRateOnset[3,:])


        indMinPvalSustain = np.unravel_index(np.argmin(pValEachCondSustain), pValEachCondSustain.shape)
        minPvalEachCellSustain[indCell] = pValEachCondSustain[indMinPvalSustain]
        minPvalIndexEachCellSustain[indCell,:] = indMinPvalSustain
        indBestSustain = np.unravel_index(np.nanargmax(np.abs(meanFiringRateSustain-firingRateEachCellBase[indCell])), meanFiringRateSustain.shape)
        bestFiringRateEachCellSustain[indCell] = meanFiringRateSustain[indBestSustain]
        bestIndexEachCellSustain[indCell,:] = indBestSustain
        maxFiringRateEachCellSustain[indCell] = np.nanmax(meanFiringRateSustain)
        minFiringRateEachCellSustain[indCell] = np.nanmin(meanFiringRateSustain)


    celldb['speechMinPvalOnset'] = minPvalEachCellOnset
    celldb['ftIndexMinPvalOnset'] = minPvalIndexEachCellOnset[:,0]
    celldb['votIndexMinPvalOnset'] = minPvalIndexEachCellOnset[:,1]
    celldb['speechMinPvalSustain'] = minPvalEachCellSustain
    celldb['ftIndexMinPvalSustain'] = minPvalIndexEachCellSustain[:,0]
    celldb['votIndexMinPvalSustain'] = minPvalIndexEachCellSustain[:,1]
    celldb['speechFiringRateBaseline'] = firingRateEachCellBase
    celldb['speechFiringRateBestOnset'] = bestFiringRateEachCellOnset
    celldb['ftIndexBestOnset'] = bestIndexEachCellOnset[:,0]
    celldb['votIndexBestOnset'] = bestIndexEachCellOnset[:,1]
    celldb['speechFiringRateBestSustain'] = bestFiringRateEachCellSustain
    celldb['ftIndexBestSustain'] = bestIndexEachCellSustain[:,0]
    celldb['votIndexBestSustain'] = bestIndexEachCellSustain[:,1]
    celldb['speechFiringRateMaxOnset'] = maxFiringRateEachCellOnset
    celldb['speechFiringRateMinOnset'] = minFiringRateEachCellOnset
    celldb['speechFiringRateMaxSustain'] = maxFiringRateEachCellSustain
    celldb['speechFiringRateMinSustain'] = minFiringRateEachCellSustain
    celldb['minFiringRate_FT_VOTmin'] = minFiringRate_FT_VOTmin
    celldb['maxFiringRate_FT_VOTmin'] = maxFiringRate_FT_VOTmin
    celldb['minFiringRate_FT_VOTmax'] = minFiringRate_FT_VOTmax
    celldb['maxFiringRate_FT_VOTmax'] = maxFiringRate_FT_VOTmax
    celldb['minFiringRate_VOT_FTmin'] = minFiringRate_VOT_FTmin
    celldb['maxFiringRate_VOT_FTmin'] = maxFiringRate_VOT_FTmin
    celldb['minFiringRate_VOT_FTmax'] = minFiringRate_VOT_FTmax
    celldb['maxFiringRate_VOT_FTmax'] = maxFiringRate_VOT_FTmax

    celldatabase.save_hdf(celldb, newdbPath)
