"""
Calculate responsiveness to each tone.
"""

import os
import sys
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
import matplotlib.pyplot as plt
from jaratoolbox import extraplots
sys.path.append('..')
import studyparams
from importlib import reload

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
fulldb = celldatabase.load_hdf(databaseFullPath)
newdbPath = os.path.join('/tmp', 'fulldb_paspeech_whichFreqResp.h5')

FIGNAME = 'selectivityIndices'
figDataFile = 'data_whichFreqs_responsive.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)

#allSubjects = studyparams.TEST_MOUSE
allSubjects = studyparams.EPHYS_MICE

N_FREQ = 16 # HARDCODED
periodsName = ['base200', 'resp100']
allPeriods = [ [-0.2, 0], [0, 0.1] ]
periodDuration = [x[1]-x[0] for x in allPeriods]

tonePvalRespEachFreq = np.full([len(fulldb), N_FREQ], np.nan)
#fulldb.tonePvalRespEachFreq = np.full([len(fulldb), N_FREQ],-1)

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse

    celldb = fulldb[fulldb.subject == thisMouse]
    nCells = len(celldb)


    meanFiringEachPeriodEachCell = np.empty((nCells, len(allPeriods)))

    #pValsEachCell = np.empty((nCells, N_FREQ))
    #minPvalEachCellResp = np.full(nCells, np.nan)
    #minPvalIndexEachCellResp = np.full(nCells, -1)
    firingRateEachCellBase = np.full(nCells, np.nan)
    #bestFiringRateEachCellResp = np.full(nCells, np.nan)
    #bestIndexEachCellResp = np.full(nCells, -1)

    #responsiveEachCell = np.empty(nCells, dtype='bool')

    indCell = -1
    for indRow, dbRow in celldb.iterrows():
        indCell += 1
        oneCell = ephyscore.Cell(dbRow)

        if 'pureTones' not in oneCell.dbRow.sessionType:
            print(f'[{indRow}] does not have a pureTones session')
            continue

        ephysData, bdata = oneCell.load('pureTones')

        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        timeRange = [-0.3, 0.6]  # In seconds

        freqEachTrial = bdata['currentFreq']
        # -- Test if trials from behavior don't match ephys --
        if (len(freqEachTrial) > len(eventOnsetTimes)) or \
           (len(freqEachTrial) < len(eventOnsetTimes)-1):
            print(f'[{indRow}] Warning! BehavTrials ({len(freqEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')
            continue
        if len(freqEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(freqEachTrial)]

        possibleFreq = np.unique(freqEachTrial)
        nFreq = len(possibleFreq)
        trialsEachCond = behavioranalysis.find_trials_each_type(freqEachTrial, possibleFreq)

        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        meanFiringEachPeriod = np.empty(len(allPeriods))
        spikesEachTrialEachPeriod = []
        for indPeriod, period in enumerate(allPeriods):
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            spikesEachTrialEachPeriod.append(spikesEachTrial)

        firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

        # NOTE: Selecting by intensity reduces the number of trials and in this case
        #       it generally yields less reliable results.
        #trialsHighIntensity = bdata['currentIntensity'] > 50
        meanFiringRateBase = np.empty(nFreq)
        meanFiringRateResp = np.empty(nFreq)
        pValEachCondResp = np.empty(nFreq)
        for indcond, thisCond in enumerate(possibleFreq):
            trialsThisCond = trialsEachCond[:,indcond]
            #trialsThisCond = trialsThisCond & trialsHighIntensity  # Select by intensity
            firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
            firingRateResp = spikesEachTrialEachPeriod[1][trialsThisCond]/periodDuration[1]
            try:
                wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateResp)
            except ValueError:
                pValThisCond = 1
            pValEachCondResp[indcond] = pValThisCond
            meanFiringRateResp[indcond] = firingRateResp.mean()

        tonePvalRespEachFreq[indRow] = np.array(pValEachCondResp)

        if indCell % 50 == 0:
            print(f'{indCell}/{nCells}')
        #print(f'[{indRow}] {str(oneCell)}')



np.savez(figDataFullPath, tonePvalRespEachFreq = tonePvalRespEachFreq)
print('saved to ' f'{figDataFullPath}')

#fulldb['tonePvalRespEachFreq'] = tonePvalRespEachFreq

#celldatabase.save_hdf(fulldb, newdbPath)
