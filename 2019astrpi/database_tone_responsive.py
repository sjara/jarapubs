"""
Calculate tone responsiveness.
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
dbPath = os.path.join(figuresDataDir, 'astrpi_laser_response.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

newdbPath = '/tmp/astrpi_tones_pval.h5'

periodsName = ['base200', 'resp100']
allPeriods = [ [-0.2, 0], [0, 0.1] ]
periodDuration = [x[1]-x[0] for x in allPeriods]
meanFiringEachPeriodEachCell = np.empty((nCells, len(allPeriods)))
N_FREQ = 16 # HARDCODED
pValsEachCell = np.empty((nCells, N_FREQ))
minPvalEachCellResp = np.full(nCells, np.nan)
minPvalIndexEachCellResp = np.full(nCells, -1)
firingRateEachCellBase = np.full(nCells, np.nan)
bestFiringRateEachCellResp = np.full(nCells, np.nan)
bestIndexEachCellResp = np.full(nCells, -1)

responsiveEachCell = np.empty(nCells, dtype='bool')

#celldb = celldb.iloc[429:]
#celldb = celldb.loc[1967:]
#celldb = celldb.loc[1995:]

indCell = -1
for indRow, dbRow in celldb.iterrows():
    #dbRow = celldb.loc[570]
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    if 'tuningCurve' not in oneCell.dbRow.sessionType:
        print(f'[{indRow}] does not have a tuningCurve session')
        continue
    
    ephysData, bdata = oneCell.load('tuningCurve')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.3, 0.6]  # In seconds

    freqEachTrial = bdata['currentFreq']
    # -- Test if trials from behavior don't match ephys --
    if (len(freqEachTrial) > len(eventOnsetTimes)) or \
       (len(freqEachTrial) < len(eventOnsetTimes)-1):
        print(f'[{indRow}] Warning! BevahTrials ({len(freqEachTrial)}) and ' +
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
    #       it generally yields higher p-values.
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
        
    indMinPvalResp = np.argmin(pValEachCondResp)
    minPvalEachCellResp[indCell] = pValEachCondResp[indMinPvalResp]
    minPvalIndexEachCellResp[indCell] = pValEachCondResp[indMinPvalResp]
    indBestResponse = np.argmax(np.abs(meanFiringRateResp-firingRateEachCellBase[indCell]))
    bestFiringRateEachCellResp[indCell] = meanFiringRateResp[indBestResponse]
    bestIndexEachCellResp[indCell] = indBestResponse
    
    if indCell % 50 == 0:
        print(f'{indCell}/{nCells}')
    #print(f'[{indRow}] {str(oneCell)}')
    
    if 0:
        correctedAlpha = 0.05/nFreq
        responsiveThisCell = np.any(pValEachCondResp<correctedAlpha)
        responsiveEachCell[indCell] = responsiveThisCell

        pValsEachCell[indCell, :] = pValEachCondResp
        pValsStr = ' '.join([f'{p:0.4f}' for p in pValEachCondResp])   
        #print(pValsStr)
        
        print(f'B:{firingRateEachCellBase[indCell]:0.2f}   ' +
              f'R:{bestFiringRateEachCellResp[indCell]:0.2f}   ' +
              f'[{bestIndexEachCellResp[indCell]}]')
        print()
        
        plt.cla()
        #timeRange = [-0.28, 0.58]
        fRaster = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                         timeRange, trialsEachCond)
        #pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset,
        #                                               indexLimitsEachTrial,timeRange)
        fontWeight = 'bold' if responsiveThisCell else None
        thisTitle = plt.title(f'[{indRow}] {minPvalEachCellResp[indCell]:0.4f}',
                              fontweight=fontWeight)
        plt.show()
        plt.waitforbuttonpress()
        #plt.pause(0.1)


celldb['toneMinPval'] = minPvalEachCellResp
celldb['toneIndexMinPval'] = minPvalIndexEachCellResp
celldb['toneFiringRateBaseline'] = firingRateEachCellBase
celldb['toneFiringRateBest'] = bestFiringRateEachCellResp
celldb['toneIndexBest'] = bestIndexEachCellResp

celldatabase.save_hdf(celldb, newdbPath)


'''
bdata.events['eventTime'][bdata.events['eventCode']==-1]
'''

'''
np.sum((celldb.toneMinPval<(0.05/16)) & (celldb.laserPvalB200R50<(0.01)))   
np.sum((celldb.toneMinPval<(0.05/16)) & (celldb.laserPvalB200R50>(0.05)))   
'''
