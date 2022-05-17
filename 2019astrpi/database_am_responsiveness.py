"""
Calculate responsiveness to amplitude modulated sounds.
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


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_tones_pval.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

newdbPath = '/tmp/astrpi_am_pval.h5'

periodsName = ['base200', 'respOnset', 'respSustained']
allPeriods = [ [-0.2, 0], [0, 0.1] , [0.1, 0.5] ]
periodDuration = [x[1]-x[0] for x in allPeriods]
meanFiringEachPeriodEachCell = np.empty((nCells, len(allPeriods)))
N_RATE = 11 # HARDCODED

pValsEachCellOnset = np.empty((nCells, N_RATE))
pValsEachCellSustain = np.empty((nCells, N_RATE))
minPvalEachCellOnset = np.full(nCells, np.nan)
minPvalEachCellSustain = np.full(nCells, np.nan)
minPvalIndexEachCellOnset = np.full(nCells, -1)
minPvalIndexEachCellSustain = np.full(nCells, -1)

firingRateEachCellBase = np.full(nCells, np.nan)
bestFiringRateEachCellOnset = np.full(nCells, np.nan)
bestFiringRateEachCellSustain = np.full(nCells, np.nan)
bestIndexEachCellOnset = np.full(nCells, -1)
bestIndexEachCellSustain = np.full(nCells, -1)
maxFiringRateEachCellOnset = np.full(nCells, np.nan)
minFiringRateEachCellOnset = np.full(nCells, np.nan)
maxFiringRateEachCellSustain = np.full(nCells, np.nan)
minFiringRateEachCellSustain = np.full(nCells, np.nan)

responsiveEachCellOnset = np.empty(nCells, dtype='bool')
responsiveEachCellSustain = np.empty(nCells, dtype='bool')

#celldb = celldb.iloc[0:]
#celldb = celldb.loc[1350:]
#celldb = celldb.loc[684:]
#celldb = celldb.loc[710:]
# Missing: 717, 723, (803), 1316

indCell = -1
for indRow, dbRow in celldb.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    if 'am' not in oneCell.dbRow.sessionType:
        print(f'[{indRow}] does not have an AM session')
        continue
    
    ephysData, bdata = oneCell.load('am')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.4, 0.8]  # In seconds

    rateEachTrial = bdata['currentFreq']
    # -- Test if trials from behavior don't match ephys --
    if (len(rateEachTrial) > len(eventOnsetTimes)) or \
       (len(rateEachTrial) < len(eventOnsetTimes)-1):
        print(f'[{indRow}] Warning! BevahTrials ({len(rateEachTrial)}) and ' +
              f'EphysTrials ({len(eventOnsetTimes)})')
        continue
    if len(rateEachTrial) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(rateEachTrial)]
    
    possibleRate = np.unique(rateEachTrial)
    nRate = len(possibleRate)
    trialsEachCond = behavioranalysis.find_trials_each_type(rateEachTrial, possibleRate)

    (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    meanFiringEachPeriod = np.empty(len(allPeriods))
    spikesEachTrialEachPeriod = []
    for indPeriod, period in enumerate(allPeriods):
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial, period)
        spikesEachTrial = spikeCountMat[:,0]
        spikesEachTrialEachPeriod.append(spikesEachTrial)

    firingRateEachCellBase[indCell] = spikesEachTrialEachPeriod[0].mean()/periodDuration[0]

    meanFiringRateBase = np.empty(nRate)
    meanFiringRateOnset = np.empty(nRate)
    pValEachCondOnset = np.empty(nRate)
    meanFiringRateSustain = np.empty(nRate)
    pValEachCondSustain = np.empty(nRate)
    for indcond, thisCond in enumerate(possibleRate):
        trialsThisCond = trialsEachCond[:,indcond]
        firingRateBase = spikesEachTrialEachPeriod[0][trialsThisCond]/periodDuration[0]
        firingRateOnset = spikesEachTrialEachPeriod[1][trialsThisCond]/periodDuration[1]
        firingRateSustain = spikesEachTrialEachPeriod[2][trialsThisCond]/periodDuration[2]
        try:
            wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateOnset)
        except ValueError:
            pValThisCond = 1
        pValEachCondOnset[indcond] = pValThisCond
        try:
            wStat, pValThisCond = stats.wilcoxon(firingRateBase, firingRateSustain)
        except ValueError:
            pValThisCond = 1
        pValEachCondSustain[indcond] = pValThisCond
        meanFiringRateOnset[indcond] = firingRateOnset.mean()
        meanFiringRateSustain[indcond] = firingRateSustain.mean()
        
    indMinPvalOnset = np.argmin(pValEachCondOnset)
    minPvalEachCellOnset[indCell] = pValEachCondOnset[indMinPvalOnset]
    minPvalIndexEachCellOnset[indCell] = indMinPvalOnset
    indBestOnset = np.argmax(np.abs(meanFiringRateOnset-firingRateEachCellBase[indCell]))
    bestFiringRateEachCellOnset[indCell] = meanFiringRateOnset[indBestOnset]
    bestIndexEachCellOnset[indCell] = indBestOnset
    maxFiringRateEachCellOnset[indCell] = np.max(meanFiringRateOnset)
    minFiringRateEachCellOnset[indCell] = np.min(meanFiringRateOnset)
    
    indMinPvalSustain = np.argmin(pValEachCondSustain)
    minPvalEachCellSustain[indCell] = pValEachCondSustain[indMinPvalSustain]
    minPvalIndexEachCellSustain[indCell] = indMinPvalSustain
    indBestSustain = np.argmax(np.abs(meanFiringRateSustain-firingRateEachCellBase[indCell]))
    bestFiringRateEachCellSustain[indCell] = meanFiringRateSustain[indBestSustain]
    bestIndexEachCellSustain[indCell] = indBestSustain
    maxFiringRateEachCellSustain[indCell] = np.max(meanFiringRateSustain)
    minFiringRateEachCellSustain[indCell] = np.min(meanFiringRateSustain)

    '''
    indMinPvalSustain = np.argmin(pValEachCondSustain)
    minPvalEachCellSustain[indCell] = pValEachCondSustain[indMinPvalSustain]
    bestFiringRateEachCellSustain[indCell] = meanFiringRateSustain[indMinPvalSustain]
    bestIndexEachCellSustain[indCell] = indMinPvalSustain
    '''
    
    if indCell % 50 == 0:
        print(f'{indCell}/{nCells}')
    #print(f'[{indRow}] {str(oneCell)}')
    
    if 0:
        correctedAlpha = 0.05/nRate
        responsiveThisCellOnset = np.any(pValEachCondOnset<correctedAlpha)
        responsiveEachCellOnset[indCell] = responsiveThisCellOnset
        responsiveThisCellSustain = np.any(pValEachCondSustain<correctedAlpha)
        responsiveEachCellSustain[indCell] = responsiveThisCellSustain


        pValsEachCellOnset[indCell, :] = pValEachCondOnset
        pValsEachCellSustain[indCell, :] = pValEachCondSustain
        pValsStr = ' '.join([f'{p:0.4f}' for p in pValEachCondOnset])   
        #print(pValsStr)
        pValsStr = ' '.join([f'{p:0.4f}' for p in pValEachCondSustain])   
        #print(pValsStr)
        print(f'B:{firingRateEachCellBase[indCell]:0.2f}   ' +
              f'O:{bestFiringRateEachCellOnset[indCell]:0.2f}   ' +
              f'S:{bestFiringRateEachCellSustain[indCell]:0.2f}')
        print()
        
        plt.cla()
        fRaster = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                         timeRange, trialsEachCond)
        #pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset,
        #                                               indexLimitsEachTrial,timeRange)
        fontWeight = 'bold' if responsiveThisCellOnset else None
        thisTitle = plt.title(f'[{indRow}] {minPvalEachCellOnset[indCell]:0.4f}   '+
                              f'{minPvalEachCellSustain[indCell]:0.4f}',
                              fontweight=fontWeight)
        plt.show()
        #plt.waitforbuttonpress()
        sys.exit()
        #plt.pause(0.1)


celldb['amMinPvalOnset'] = minPvalEachCellOnset
celldb['amIndexMinPvalOnset'] = minPvalIndexEachCellOnset
celldb['amMinPvalSustain'] = minPvalEachCellSustain
celldb['amIndexMinPvalOnset'] = minPvalIndexEachCellSustain

celldb['amFiringRateBaseline'] = firingRateEachCellBase
celldb['amFiringRateBestOnset'] = bestFiringRateEachCellOnset
celldb['amIndexBestOnset'] = bestIndexEachCellOnset
celldb['amFiringRateBestSustain'] = bestFiringRateEachCellSustain
celldb['amIndexBestSustain'] = bestIndexEachCellSustain

celldb['amFiringRateMaxOnset'] = maxFiringRateEachCellOnset
celldb['amFiringRateMinOnset'] = minFiringRateEachCellOnset
celldb['amFiringRateMaxSustain'] = maxFiringRateEachCellSustain
celldb['amFiringRateMinSustain'] = minFiringRateEachCellSustain

celldatabase.save_hdf(celldb, newdbPath)


'''
np.sum((celldb.amMinPval<(0.05/11)) & (celldb.laserPvalB200R50<(0.01)))   
np.sum((celldb.amMinPval<(0.05/11)) & (celldb.laserPvalB200R50>(0.05)))   
'''

