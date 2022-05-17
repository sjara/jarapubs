"""
Calculate laser response to see if cells are D1 or nD1.
"""

import os
import sys
import studyparams
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from scipy import stats

import matplotlib.pyplot as plt
from jaratoolbox import extraplots


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_included_cells.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

newdbPath = '/tmp/astrpi_laser_response.h5'

#periodsName = ['base50', 'resp50', 'base100', 'resp100', 'base200']
#allPeriods = [ [-0.05, 0], [0, 0.05], [-0.1, 0], [0, 0.1] , [-0.2, 0] ]
periodsName = ['base200', 'resp50']
allPeriods = [ [-0.2, 0], [0, 0.05] ]
periodDuration = [x[1]-x[0] for x in allPeriods]
meanFiringEachPeriodEachCell = np.empty((nCells, len(allPeriods)))
#pValEachCell100 = np.empty(nCells)
#pValEachCell50 = np.empty(nCells)
#pValEachCellB100R50 = np.empty(nCells)
pValEachCellB200R50 = np.empty(nCells)

indCell = -1
for indRow, dbRow in celldb.iterrows():
    indCell += 1
    oneCell = ephyscore.Cell(dbRow)

    ephysData, bdata = oneCell.load('laserpulse')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.3, 0.6]  # In seconds

    (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    meanFiringEachPeriod = np.empty(len(allPeriods))
    spikesEachTrialEachPeriod = []
    for indPeriod, period in enumerate(allPeriods):
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial, period)
        spikesEachTrial = spikeCountMat[:,0]
        spikesEachTrialEachPeriod.append(spikesEachTrial)
        
    firingRateEachTrialEachPeriod = np.array(spikesEachTrialEachPeriod) / \
                                    np.array(periodDuration)[:,np.newaxis]
    meanFiringEachPeriodEachCell[indCell,:] = firingRateEachTrialEachPeriod.mean(axis=1)
    
    try:
        wStat, pValB200R50 = stats.wilcoxon(firingRateEachTrialEachPeriod[0,:],
                                            firingRateEachTrialEachPeriod[1,:])
    except ValueError:
        pValB200R50 = 1

    if indCell % 100 == 0:
        print(f'{indCell}/{nCells}')
    #print(f'[{indRow}] {str(oneCell)}')
    pValEachCellB200R50[indCell] = pValB200R50
    
    if 0:
        plt.cla()
        pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset,
                                                       indexLimitsEachTrial,timeRange)
        plt.title(f'[{indRow}] p = {pValB200R50:0.4f}   ' +
                  f'{meanFiringEachPeriodEachCell[indCell,0]:0.3f} ' +
                  f'vs {meanFiringEachPeriodEachCell[indCell,1]:0.3f}')
        plt.show()
        plt.waitforbuttonpress()
        #plt.pause(0.1)

celldb['laserPvalB200R50'] = pValEachCellB200R50
celldb['laserBaseRate200'] = meanFiringEachPeriodEachCell[:,0]
celldb['laserRespRate50'] = meanFiringEachPeriodEachCell[:,1]

celldatabase.save_hdf(celldb, newdbPath)

