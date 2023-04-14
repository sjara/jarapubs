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

FIGNAME = 'selectivityIndices'
figDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
allSubjects = studyparams.EPHYS_MICE
#allSubjects = studyparams.TEST_MOUSE

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
nCells = len(celldb)

#shuffledSI = np.empty([nCells, 2000])
shuffledSIVot_FTmin = np.empty([nCells, 2000])
shuffledSIVot_FTmax = np.empty([nCells, 2000])
shuffledSIFt_VOTmin = np.empty([nCells, 2000])
shuffledSIFt_VOTmax = np.empty([nCells, 2000])
#shuffledMinRate = np.empty([nCells, 2000])
#shuffledMaxRate = np.empty([nCells, 2000])
shuffledMinRateVot_FTmin = np.empty([nCells, 2000])
shuffledMaxRateVot_FTmin = np.empty([nCells, 2000])
shuffledMinRateVot_FTmax = np.empty([nCells, 2000])
shuffledMaxRateVot_FTmax = np.empty([nCells, 2000])
shuffledMinRateFt_VOTmin = np.empty([nCells, 2000])
shuffledMaxRateFt_VOTmin = np.empty([nCells, 2000])
shuffledMinRateFt_VOTmax = np.empty([nCells, 2000])
shuffledMaxRateFt_VOTmax = np.empty([nCells, 2000])
#votSI = np.empty(nCells)
#minRate = np.empty(nCells)
#maxRate = np.empty(nCells)
minRateVotFTmin = np.empty(nCells)
maxRateVotFTmin = np.empty(nCells)
votSI_FTmin = np.empty(nCells)
minRateVotFTmax = np.empty(nCells)
maxRateVotFTmax = np.empty(nCells)
votSI_FTmax = np.empty(nCells)
minRateFtVOTmin = np.empty(nCells)
maxRateFtVOTmin = np.empty(nCells)
ftSI_VOTmin = np.empty(nCells)
minRateFtVOTmax = np.empty(nCells)
maxRateFtVOTmax = np.empty(nCells)
ftSI_VOTmax = np.empty(nCells)
pvalPermutationtestVot = np.ones(nCells)
pvalPermutationtestFt = np.ones(nCells)
alpha = 0.025

indCell = -1
for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    #dbPath = os.path.join(databaseDir, f'{subject}_paspeech_am_tuning.h5')
    #celldb = celldatabase.load_hdf(dbPath)
    mousedb = celldb[celldb.subject==thisMouse]
    nCellsThisMouse = len(mousedb)

    #newdbPath = os.path.join(databaseDir, f'{subject}t.h5')
    '''
    periodsName = ['base200', 'respOnset', 'respSustained']
    allPeriods = [ [-0.2, 0], [0, 0.12] , [0.12, 0.24] ] #try with shorter period for onset response.
    '''
    periodsName = ['respOnset']
    allPeriods = [[0, 0.12]]
    #for speed and usefulness, only calculating shuffled indices for response onset (which is what we used for the selectivity indices)
    periodDuration = [x[1] - x[0] for x in allPeriods]

    for indRow, dbRow in mousedb.iterrows():
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

        ## -- Test if trials from behavior don't match ephys --
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

        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        spikesEachTrialEachPeriod = []
        for indPeriod, period in enumerate(allPeriods):
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            spikesEachTrialEachPeriod.append(spikesEachTrial)

            ## Calculate mean firing rates and responsiveness for each speech sound (FT-VOT combination)
            #-- Calculate actual SI *NOT collapsing across irrelevant feature*
            trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)

            meanFiringRateOnset = np.empty([nFT, nVOT])

            for indFT, thisFT in enumerate(possibleFTParams):
                for indVOT, thisVOT in enumerate(possibleVOTParams):
                    trialsThisCond = trialsEachCond[:, indFT, indVOT]
                    firingRateOnset = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
                    meanFiringRateOnset[indFT, indVOT] = firingRateOnset.mean()

            minFiringRateVot_FTmin = np.min(meanFiringRateOnset[0,:])
            maxFiringRateVot_FTmin = np.max(meanFiringRateOnset[0,:])
            minRateVotFTmin[indCell] = minFiringRateVot_FTmin
            maxRateVotFTmin[indCell] = maxFiringRateVot_FTmin
            votSI_FTmin[indCell] = (maxFiringRateVot_FTmin - minFiringRateVot_FTmin)/(maxFiringRateVot_FTmin + minFiringRateVot_FTmin)

            minFiringRateVot_FTmax = np.min(meanFiringRateOnset[3,:])
            maxFiringRateVot_FTmax = np.max(meanFiringRateOnset[3,:])
            minRateVotFTmax[indCell] = minFiringRateVot_FTmax
            maxRateVotFTmax[indCell] = maxFiringRateVot_FTmax
            votSI_FTmax[indCell] = (maxFiringRateVot_FTmax - minFiringRateVot_FTmax)/(maxFiringRateVot_FTmax + minFiringRateVot_FTmax)

            minFiringRateFt_VOTmin = np.min(meanFiringRateOnset[:,0])
            maxFiringRateFt_VOTmin = np.max(meanFiringRateOnset[:,0])
            minRateFtVOTmin[indCell] = minFiringRateFt_VOTmin
            maxRateFtVOTmin[indCell] = maxFiringRateFt_VOTmin
            ftSI_VOTmin[indCell] = (maxFiringRateFt_VOTmin - minFiringRateFt_VOTmin)/(maxFiringRateFt_VOTmin + minFiringRateFt_VOTmin)

            minFiringRateFt_VOTmax = np.min(meanFiringRateOnset[:,3])
            maxFiringRateFt_VOTmax = np.max(meanFiringRateOnset[:,3])
            minRateFtVOTmax[indCell] = minFiringRateFt_VOTmax
            maxRateFtVOTmax[indCell] = maxFiringRateFt_VOTmax
            ftSI_VOTmax[indCell] = (maxFiringRateFt_VOTmax - minFiringRateFt_VOTmax)/(maxFiringRateFt_VOTmax + minFiringRateFt_VOTmax)


            trialsEachVot = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)
            shuffledVOTParamsEachTrial = VOTParamsEachTrial.copy()
            trialsEachFt = behavioranalysis.find_trials_each_type(FTParamsEachTrial, possibleFTParams)
            shuffledFTParamsEachTrial = FTParamsEachTrial.copy()

            for indShuffle, thisShuffle in enumerate(shuffledSI[0]):

                np.random.shuffle(shuffledVOTParamsEachTrial)
                trialsEachCond_VOTshuffled = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, shuffledVOTParamsEachTrial, possibleVOTParams)
                #--FT
                trialsEachCond_FTshuffled = behavioranalysis.find_trials_each_combination(shuffledFTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)

                meanFiringRateOnsetShuffleVot = np.empty([nFT, nVOT])
                meanFiringRateOnsetShuffleFt = np.empty([nFT, nVOT])

                for indFT, thisFT in enumerate(possibleFTParams):
                    for indVOT, thisVOT in enumerate(possibleVOTParams):
                        #trialsThisCond = trialsEachCond[:, indFT, indVOT]
                        trialsThisVotShuffled = trialsEachCond_VOTshuffled[:,indFT,indVOT]
                        trialsThisFtShuffled = trialsEachCond_FTshuffled[:,indFT,indVOT]

                        firingRateOnsetVotShuffle = spikesEachTrialEachPeriod[0][trialsThisVotShuffled]/periodDuration
                        meanFiringRateOnsetShuffleVot[indFT, indVOT] = firingRateOnsetVotShuffle.mean()

                        firingRateOnsetFtShuffle = spikesEachTrialEachPeriod[0][trialsThisFtShuffled]/periodDuration
                        meanFiringRateOnsetShuffleFt[indFT, indVOT] = firingRateOnsetFtShuffle.mean()

                minFiringRateVot_FTmin = np.min(meanFiringRateOnsetShuffleVot[0,:])
                maxFiringRateVot_FTmin = np.max(meanFiringRateOnsetShuffleVot[0,:])
                shuffledMinRateVot_FTmin[indCell, indShuffle] = minFiringRateVot_FTmin
                shuffledMaxRateVot_FTmin[indCell, indShuffle] = maxFiringRateVot_FTmin
                shuffledSIVot_FTmin[indCell, indShuffle] = (maxFiringRateVot_FTmin - minFiringRateVot_FTmin)/(maxFiringRateVot_FTmin + minFiringRateVot_FTmin)

                minFiringRateVot_FTmax = np.min(meanFiringRateOnsetShuffleVot[-1,:])
                maxFiringRateVot_FTmax = np.max(meanFiringRateOnsetShuffleVot[-1,:])
                shuffledMinRateVot_FTmax[indCell, indShuffle] = minFiringRateVot_FTmax
                shuffledMaxRateVot_FTmax[indCell, indShuffle] = maxFiringRateVot_FTmax
                shuffledSIVot_FTmax[indCell, indShuffle] = (maxFiringRateVot_FTmax - minFiringRateVot_FTmax)/(maxFiringRateVot_FTmax + minFiringRateVot_FTmax)


                minFiringRateFt_VOTmin = np.min(meanFiringRateOnsetShuffleFt[:,0])
                maxFiringRateFt_VOTmin = np.max(meanFiringRateOnsetShuffleFt[:,0])
                shuffledMinRateFt_VOTmin[indCell, indShuffle] = minFiringRateFt_VOTmin
                shuffledMaxRateFt_VOTmin[indCell, indShuffle] = maxFiringRateFt_VOTmin
                shuffledSIFt_VOTmin[indCell, indShuffle] = (maxFiringRateFt_VOTmin - minFiringRateFt_VOTmin)/(maxFiringRateFt_VOTmin + minFiringRateFt_VOTmin)

                minFiringRateFt_VOTmax = np.min(meanFiringRateOnsetShuffleFt[:,-1])
                maxFiringRateFt_VOTmax = np.max(meanFiringRateOnsetShuffleFt[:,-1])
                shuffledMinRateFt_VOTmax[indCell, indShuffle] = minFiringRateFt_VOTmax
                shuffledMaxRateFt_VOTmax[indCell, indShuffle] = maxFiringRateFt_VOTmax
                shuffledSIFt_VOTmax[indCell, indShuffle] = (maxFiringRateFt_VOTmax - minFiringRateFt_VOTmax)/(maxFiringRateFt_VOTmax + minFiringRateFt_VOTmax)
                '''
                if votSI_FTmax[indCell] > votSI_FTmin[indCell]:
                    pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmax[indCell,:] >=votSI_FTmax[indCell])
                elif votSI_FTmax[indCell] < votSI_FTmin[indCell]:
                    pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmin[indCell,:] >=votSI_FTmin[indCell])

                if ftSI_VOTmax[indCell] > ftSI_VOTmin[indCell]:
                    pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmax[indCell,:] >=ftSI_VOTmax[indCell])
                elif ftSI_VOTmax[indCell] < ftSI_VOTmin[indCell]:
                    pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmin[indCell,:] >=ftSI_VOTmin[indCell])
                '''

print(f'n cells w/pvalPermutationtestVot < 0.025 = {np.sum(pvalPermutationtestVot < 0.025)}')
print(f'n cells w/pvalPermutationtestFt < 0.025 = {np.sum(pvalPermutationtestFt < 0.025)}')
#print(f'n cells w/pvalPermutationtestVot < 0.05 = {np.sum(pvalPermutationtestVot < 0.05)}')
#print(f'n cells w/pvalPermutationtestFt < 0.05 = {np.sum(pvalPermutationtestFt < 0.05)}')


np.savez(figDataFullPath, shuffledSIVot_FTmin = shuffledSIVot_FTmin, shuffledSIVot_FTmax = shuffledSIVot_FTmax, shuffledSIFt_VOTmin = shuffledSIFt_VOTmin, shuffledSIFt_VOTmax = shuffledSIFt_VOTmax, shuffledMinRateVot_FTmin = shuffledMinRateVot_FTmin, shuffledMaxRateVot_FTmin = shuffledMaxRateVot_FTmin, shuffledMinRateVot_FTmax = shuffledMinRateVot_FTmax, shuffledMaxRateVot_FTmax = shuffledMaxRateVot_FTmax, shuffledMinRateFt_VOTmin = shuffledMinRateFt_VOTmin, shuffledMaxRateFt_VOTmin = shuffledMaxRateFt_VOTmin, shuffledMinRateFt_VOTmax = shuffledMinRateFt_VOTmax, shuffledMaxRateFt_VOTmax = shuffledMaxRateFt_VOTmax, minRateVotFTmin = minRateVotFTmin,  maxRateVotFTmin = maxRateVotFTmin, votSI_FTmin = votSI_FTmin, minRateVotFTmax = minRateVotFTmax, maxRateVotFTmax = maxRateVotFTmax, votSI_FTmax = votSI_FTmax, minRateFtVOTmin = minRateFtVOTmin, maxRateFtVOTmin = maxRateFtVOTmin, ftSI_VOTmin = ftSI_VOTmin, minRateFtVOTmax = minRateFtVOTmax, maxRateFtVOTmax = maxRateFtVOTmax, ftSI_VOTmax = ftSI_VOTmax, pvalPermutationtestVot = pvalPermutationtestVot, pvalPermutationtestFt = pvalPermutationtestFt, alpha = alpha)
print('saved to ' f'{figDataFullPath}')
