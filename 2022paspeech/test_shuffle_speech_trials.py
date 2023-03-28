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
#allSubjects = ['feat009', 'feat010']
#allSubjects = studyparams.TEST_MOUSE

shuffledSI = []

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    dbPath = os.path.join(databaseDir, f'{subject}_paspeech_am_tuning.h5')
    celldb = celldatabase.load_hdf(dbPath)
    nCells = len(celldb)

    #newdbPath = os.path.join(databaseDir, f'{subject}t.h5')
    '''
    periodsName = ['base200', 'respOnset', 'respSustained']
    allPeriods = [ [-0.2, 0], [0, 0.12] , [0.12, 0.24] ] #try with shorter period for onset response.
    '''
    periodsName = ['respOnset']
    allPeriods = [[0, 0.12]]
    #for speed and usefulness, only calculating shuffled indices for response onset (which is what we used for the selectivity indices)

    periodDuration = [x[1] - x[0] for x in allPeriods]


    shuffledSI_thisMouse = np.zeros([nCells, 2000])
    firingRateEachCellBase = np.full(nCells, np.nan)
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

        meanFiringEachPeriod = np.empty(len(allPeriods))
        spikesEachTrialEachPeriod = []
        for indPeriod, period in enumerate(allPeriods):
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            spikesEachTrialEachPeriod.append(spikesEachTrial)

            firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

            ## Calculate mean firing rates and responsiveness for each speech sound (FT-VOT combination)

            #meanFiringRateBase = np.empty([nFT, nVOT])
            #meanFiringRateOnset = np.empty([nFT, nVOT])
            #pValEachCondOnset = np.ones([nFT, nVOT])
            #meanFiringRateSustain = np.empty([nFT, nVOT])
            #pValEachCondSustain = np.ones([nFT, nVOT])
            trialsEachCond = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)
            shuffledVOTParamsEachTrial = VOTParamsEachTrial.copy()
            for indShuffle, thisShuffle in enumerate(shuffledSI_thisMouse[0]):
                meanFiringRateOnset = np.empty(nVOT)
                np.random.shuffle(shuffledVOTParamsEachTrial)
                shuffledTrialsEachCond = behavioranalysis.find_trials_each_type(shuffledVOTParamsEachTrial, possibleVOTParams)

                for indVOT, thisVOT in enumerate(possibleVOTParams):
                    shuffledTrialsThisCond = shuffledTrialsEachCond[:, indVOT]
                    #firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
                    firingRateOnset = spikesEachTrialEachPeriod[0][shuffledTrialsThisCond]/periodDuration
                    #firingRateSustain = spikesEachTrialEachPeriod[2][trialsThisCond]/periodDuration[2]

                    meanFiringRateOnset[indVOT] = firingRateOnset.mean()


                shuffledSI_thisMouse[indCell, indShuffle] = (np.nanmax(meanFiringRateOnset) - np.nanmin(meanFiringRateOnset))/(np.nanmax(meanFiringRateOnset) + np.nanmin(meanFiringRateOnset))

    shuffledSI.append(shuffledSI_thisMouse)
