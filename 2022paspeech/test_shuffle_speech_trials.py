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

shuffledSI = []
shuffledSI_Vot_FTmin = []
shuffledSI_Vot_FTmax = []
shuffledSI_Ft_VOTmin = []
shuffledSI_Ft_VOTmax = []

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
nCells = len(celldb)

shuffledSI = np.empty([nCells, 2000])
shuffledSIVot_FTmin = np.empty([nCells, 2000])
shuffledSIVot_FTmax = np.empty([nCells, 2000])
shuffledSIFt_VOTmin = np.empty([nCells, 2000])
shuffledSIFt_VOTmax = np.empty([nCells, 2000])
shuffledMinRate = np.empty([nCells, 2000])
shuffledMaxRate = np.empty([nCells, 2000])
shuffledMinRateVot_FTmin = np.empty([nCells, 2000])
shuffledMaxRateVot_FTmin = np.empty([nCells, 2000])
shuffledMinRateVot_FTmax = np.empty([nCells, 2000])
shuffledMaxRateVot_FTmax = np.empty([nCells, 2000])
shuffledMinRateFt_VOTmin = np.empty([nCells, 2000])
shuffledMaxRateFt_VOTmin = np.empty([nCells, 2000])
shuffledMinRateFt_VOTmax = np.empty([nCells, 2000])
shuffledMaxRateFt_VOTmax = np.empty([nCells, 2000])
votSI = np.empty(nCells)
minRate = np.empty(nCells)
maxRate = np.empty(nCells)
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

            #firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

            ## Calculate mean firing rates and responsiveness for each speech sound (FT-VOT combination)
            '''
            trialsEachVot = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)

            shuffledVOTParamsEachTrial = VOTParamsEachTrial.copy()
            meanFiringRateOnset = np.empty(nVOT)
            #-- Calculate actual VOT SI for cell *collapsing across all FTs*
            for indVOT, thisVOT in enumerate(possibleVOTParams):
                trialsThisVOT = trialsEachVot[:, indVOT]
                firingRateOnset = spikesEachTrialEachPeriod[0][trialsThisVOT]/periodDuration
                meanFiringRateOnset[indVOT] = firingRateOnset.mean()
            minFiringRateVot = np.min(meanFiringRateOnset)
            maxFiringRateVot = np.max(meanFiringRateOnset)
            minRate[indCell] = minFiringRateVot
            maxRate[indCell] = maxFiringRateVot
            votSI[indCell] = (maxFiringRateVot - minFiringRateVot)/(maxFiringRateVot + minFiringRateVot)
            '''
            #-- Calculate actual SI *NOT collapsing across irrelevant feature*
            trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)
            meanFiringRateOnsetVot_FTmin = np.empty(nVOT)
            meanFiringRateOnsetVot_FTmax = np.empty(nVOT)
            meanFiringRateOnsetFt_VOTmin = np.empty(nFT)
            meanFiringRateOnsetFt_VOTmax = np.empty(nFT)
            for indCond, thisCond in enumerate(possibleVOTParams):
                trialsThisCondVot_FTmin = trialsEachCond[:,0,indCond]
                trialsThisCondVot_FTmax = trialsEachCond[:,3,indCond]
                firingRateOnsetVot_FTmin = spikesEachTrialEachPeriod[0][trialsThisCondVot_FTmin]/periodDuration
                firingRateOnsetVot_FTmax = spikesEachTrialEachPeriod[0][trialsThisCondVot_FTmax]/periodDuration
                meanFiringRateOnsetVot_FTmin[indCond] = firingRateOnsetVot_FTmin.mean()
                meanFiringRateOnsetVot_FTmax[indCond] = firingRateOnsetVot_FTmax.mean()

                trialsThisCondFt_VOTmin = trialsEachCond[:,indCond,0]
                trialsThisCondFt_VOTmax = trialsEachCond[:,indCond,3]
                firingRateOnsetFt_VOTmin = spikesEachTrialEachPeriod[0][trialsThisCondFt_VOTmin]/periodDuration
                firingRateOnsetFt_VOTmax = spikesEachTrialEachPeriod[0][trialsThisCondFt_VOTmax]/periodDuration
                meanFiringRateOnsetVot_FTmin[indCond] = firingRateOnsetFt_VOTmin.mean()
                meanFiringRateOnsetVot_FTmax[indCond] = firingRateOnsetFt_VOTmax.mean()

            minFiringRateVot_FTmin = np.min(meanFiringRateOnsetVot_FTmin)
            maxFiringRateVot_FTmin = np.max(meanFiringRateOnsetVot_FTmin)
            minRateVotFTmin[indCell] = minFiringRateVot_FTmin
            maxRateVotFTmin[indCell] = maxFiringRateVot_FTmin
            votSI_FTmin[indCell] = (maxFiringRateVot_FTmin - minFiringRateVot_FTmin)/(maxFiringRateVot_FTmin + minFiringRateVot_FTmin)

            minFiringRateVot_FTmax = np.min(meanFiringRateOnsetVot_FTmax)
            maxFiringRateVot_FTmax = np.max(meanFiringRateOnsetVot_FTmax)
            minRateVotFTmax[indCell] = minFiringRateVot_FTmax
            maxRateVotFTmax[indCell] = maxFiringRateVot_FTmax
            votSI_FTmax[indCell] = (maxFiringRateVot_FTmax - minFiringRateVot_FTmax)/(maxFiringRateVot_FTmax + minFiringRateVot_FTmax)

            minFiringRateFt_VOTmin = np.min(meanFiringRateOnsetFt_VOTmin)
            maxFiringRateFt_VOTmin = np.max(meanFiringRateOnsetFt_VOTmin)
            minRateFtVOTmin[indCell] = minFiringRateFt_VOTmin
            maxRateFtVOTmin[indCell] = maxFiringRateFt_VOTmin
            ftSI_VOTmin[indCell] = (maxFiringRateFt_VOTmin - minFiringRateFt_VOTmin)/(maxFiringRateFt_VOTmin + minFiringRateFt_VOTmin)

            minFiringRateFt_VOTmax = np.min(meanFiringRateOnsetFt_VOTmax)
            maxFiringRateFt_VOTmax = np.max(meanFiringRateOnsetFt_VOTmax)
            minRateFtVOTmax[indCell] = minFiringRateFt_VOTmax
            maxRateFtVOTmax[indCell] = maxFiringRateFt_VOTmax
            ftSI_VOTmax[indCell] = (maxFiringRateFt_VOTmax - minFiringRateFt_VOTmax)/(maxFiringRateFt_VOTmax + minFiringRateFt_VOTmax)

            trialsEachVot = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)
            shuffledVOTParamsEachTrial = VOTParamsEachTrial.copy()
            trialsEachFt = behavioranalysis.find_trials_each_type(FTParamsEachTrial, possibleFTParams)
            shuffledFTParamsEachTrial = FTParamsEachTrial.copy()
            for indShuffle, thisShuffle in enumerate(shuffledSI[0]):
                #meanFiringRateOnset = np.empty(nVOT)
                meanFiringRateOnsetVot_FTmin = np.empty(nVOT)
                meanFiringRateOnsetVot_FTmax = np.empty(nVOT)
                meanFiringRateOnsetVot_FTmin = np.empty(nVOT)
                meanFiringRateOnsetVot_FTmax = np.empty(nVOT)
                np.random.shuffle(shuffledVOTParamsEachTrial)
                shuffledtrialsEachVot = behavioranalysis.find_trials_each_type(shuffledVOTParamsEachTrial, possibleVOTParams)
                trialsEachCond_VOTshuffled = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, shuffledVOTParamsEachTrial, possibleVOTParams)
                #--FT
                meanFiringRateOnsetFt_VOTmin = np.empty(nFT)
                meanFiringRateOnsetFt_VOTmax = np.empty(nFT)
                meanFiringRateOnsetFt_VOTmin = np.empty(nFT)
                meanFiringRateOnsetFt_VOTmax = np.empty(nFT)
                np.random.shuffle(shuffledFTParamsEachTrial)
                shuffledtrialsEachFt = behavioranalysis.find_trials_each_type(shuffledFTParamsEachTrial, possibleFTParams)
                trialsEachCond_FTshuffled = behavioranalysis.find_trials_each_combination(shuffledFTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)

                for indVOT, thisVOT in enumerate(possibleVOTParams):
                    #shuffledtrialsThisVOT = shuffledtrialsEachVot[:, indVOT]
                    #firingRateOnset = spikesEachTrialEachPeriod[0][shuffledtrialsThisVOT]/periodDuration
                    #meanFiringRateOnset[indVOT] = firingRateOnset.mean()
                    shuffledtrialsThisVOT_FTmin = trialsEachCond_VOTshuffled[:,0,indVOT]
                    shuffledtrialsThisVOT_FTmax = trialsEachCond_VOTshuffled[:,3,indVOT]
                    firingRateOnsetVot_FTmin = spikesEachTrialEachPeriod[0][shuffledtrialsThisVOT_FTmin]/periodDuration
                    meanFiringRateOnsetVot_FTmin[indVOT] = firingRateOnsetVot_FTmin.mean()
                    firingRateOnsetVot_FTmax = spikesEachTrialEachPeriod[0][shuffledtrialsThisVOT_FTmax]/periodDuration
                    meanFiringRateOnsetVot_FTmax[indVOT] = firingRateOnsetVot_FTmax.mean()

                    shuffledtrialsThisFT_VOTmin = trialsEachCond_FTshuffled[:,indVOT,0]
                    shuffledtrialsThisFT_VOTmax = trialsEachCond_FTshuffled[:,indVOT,3]
                    firingRateOnsetFt_VOTmin = spikesEachTrialEachPeriod[0][shuffledtrialsThisFT_VOTmin]/periodDuration
                    meanFiringRateOnsetFt_VOTmin[indVOT] = firingRateOnsetFt_VOTmin.mean()
                    firingRateOnsetFt_VOTmax = spikesEachTrialEachPeriod[0][shuffledtrialsThisFT_VOTmax]/periodDuration
                    meanFiringRateOnsetFt_VOTmax[indVOT] = firingRateOnsetFt_VOTmax.mean()
                '''
                minFiringRateVot = np.min(meanFiringRateOnset)
                maxFiringRateVot = np.max(meanFiringRateOnset)
                shuffledMinRate[indCell, indShuffle] = minFiringRateVot
                shuffledMaxRate[indCell, indShuffle] = maxFiringRateVot

                shuffledSI[indCell, indShuffle] = (maxFiringRateVot - minFiringRateVot)/(maxFiringRateVot + minFiringRateVot)
                '''
                minFiringRateVot_FTmin = np.min(meanFiringRateOnsetVot_FTmin)
                maxFiringRateVot_FTmin = np.max(meanFiringRateOnsetVot_FTmin)
                shuffledMinRateVot_FTmin[indCell, indShuffle] = minFiringRateVot_FTmin
                shuffledMaxRateVot_FTmin[indCell, indShuffle] = maxFiringRateVot_FTmin
                shuffledSIVot_FTmin[indCell, indShuffle] = (maxFiringRateVot_FTmin - minFiringRateVot_FTmin)/(maxFiringRateVot_FTmin + minFiringRateVot_FTmin)

                minFiringRateVot_FTmax = np.min(meanFiringRateOnsetVot_FTmax)
                maxFiringRateVot_FTmax = np.max(meanFiringRateOnsetVot_FTmax)
                shuffledMinRateVot_FTmax[indCell, indShuffle] = minFiringRateVot_FTmax
                shuffledMaxRateVot_FTmax[indCell, indShuffle] = maxFiringRateVot_FTmax
                shuffledSIVot_FTmax[indCell, indShuffle] = (maxFiringRateVot_FTmax - minFiringRateVot_FTmax)/(maxFiringRateVot_FTmax + minFiringRateVot_FTmax)


                minFiringRateFt_VOTmin = np.min(meanFiringRateOnsetFt_VOTmin)
                maxFiringRateFt_VOTmin = np.max(meanFiringRateOnsetFt_VOTmin)
                shuffledMinRateFt_VOTmin[indCell, indShuffle] = minFiringRateFt_VOTmin
                shuffledMaxRateFt_VOTmin[indCell, indShuffle] = maxFiringRateFt_VOTmin
                shuffledSIFt_VOTmin[indCell, indShuffle] = (maxFiringRateFt_VOTmin - minFiringRateFt_VOTmin)/(maxFiringRateFt_VOTmin + minFiringRateFt_VOTmin)

                minFiringRateFt_VOTmax = np.min(meanFiringRateOnsetFt_VOTmax)
                maxFiringRateFt_VOTmax = np.max(meanFiringRateOnsetFt_VOTmax)
                shuffledMinRateFt_VOTmax[indCell, indShuffle] = minFiringRateFt_VOTmax
                shuffledMaxRateFt_VOTmax[indCell, indShuffle] = maxFiringRateFt_VOTmax
                shuffledSIFt_VOTmax[indCell, indShuffle] = (maxFiringRateFt_VOTmax - minFiringRateFt_VOTmax)/(maxFiringRateFt_VOTmax + minFiringRateFt_VOTmax)

                if votSI_FTmax[indCell] > votSI_FTmin[indCell]:
                    pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmax[indCell,:] >=votSI_FTmax[indCell])
                elif votSI_FTmax[indCell] < votSI_FTmin[indCell]:
                    pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmin[indCell,:] >=votSI_FTmin[indCell])

                if ftSI_VOTmax[indCell] > ftSI_VOTmin[indCell]:
                    pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmax[indCell,:] >=ftSI_VOTmax[indCell])
                elif ftSI_VOTmax[indCell] < ftSI_VOTmin[indCell]:
                    pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmin[indCell,:] >=ftSI_VOTmin[indCell])

print(f'n cells w/pvalPermutationtestVot < 0.025 = {np.sum(pvalPermutationtestVot < 0.025)}')
print(f'n cells w/pvalPermutationtestFt < 0.025 = {np.sum(pvalPermutationtestFt < 0.025)}')
#print(f'n cells w/pvalPermutationtestVot < 0.05 = {np.sum(pvalPermutationtestVot < 0.05)}')
#print(f'n cells w/pvalPermutationtestFt < 0.05 = {np.sum(pvalPermutationtestFt < 0.05)}')


np.savez(figDataFullPath, shuffledSIVot_FTmin = shuffledSIVot_FTmin, shuffledSIVot_FTmax = shuffledSIVot_FTmax, shuffledSIFt_VOTmin = shuffledSIFt_VOTmin, shuffledSIFt_VOTmax = shuffledSIFt_VOTmax, shuffledMinRateVot_FTmin = shuffledMinRateVot_FTmin, shuffledMaxRateVot_FTmin = shuffledMaxRateVot_FTmin, shuffledMinRateVot_FTmax = shuffledMinRateVot_FTmax, shuffledMaxRateVot_FTmax = shuffledMaxRateVot_FTmax, shuffledMinRateFt_VOTmin = shuffledMinRateFt_VOTmin, shuffledMaxRateFt_VOTmin = shuffledMaxRateFt_VOTmin, shuffledMinRateFt_VOTmax = shuffledMinRateFt_VOTmax, shuffledMaxRateFt_VOTmax = shuffledMaxRateFt_VOTmax, minRateVotFTmin = minRateVotFTmin,  maxRateVotFTmin = maxRateVotFTmin, votSI_FTmin = votSI_FTmin, minRateVotFTmax = minRateVotFTmax, maxRateVotFTmax = maxRateVotFTmax, votSI_FTmax = votSI_FTmax, minRateFtVOTmin = minRateFtVOTmin, maxRateFtVOTmin = maxRateFtVOTmin, ftSI_VOTmin = ftSI_VOTmin, minRateFtVOTmax = minRateFtVOTmax, maxRateFtVOTmax = maxRateFtVOTmax, ftSI_VOTmax = ftSI_VOTmax, pvalPermutationtestVot = pvalPermutationtestVot, pvalPermutationtestFt = pvalPermutationtestFt, alpha = alpha)
print('saved to ' f'{figDataFullPath}')
