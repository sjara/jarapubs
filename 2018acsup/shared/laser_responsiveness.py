"""
Update the database with columns related to responsiveness to laser.

This takes about 1.5 min.

Based on database_generation_funcs.laser_response()
"""

import numpy as np
from scipy import stats
from jaratoolbox import celldatabase
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
import matplotlib.pyplot as plt

timeRange = [-0.5, 1]  # In seconds

# -- For response in first 10 ms --
baseRange10ms = [-0.05, -0.04]
responseRange10ms = [0.0, 0.01]

# -- For response in first 50 ms --
baseRange50ms = [-0.1, -0.05]
responseRange50ms = [0.0, 0.05]

dbPath = '/data/figuresdata/2018acsup/shared/celldb_2018acsup_photoid_basic.h5'
dbPathNew = '/data/figuresdata/2018acsup/shared/celldb_2018acsup_photoid.h5'
celldb = celldatabase.load_hdf(dbPath)

baseFiringRateAll10ms = np.empty(len(celldb))
laserFiringRateAll10ms = np.empty(len(celldb))
pValLaserAll10ms = np.empty(len(celldb))
baseFiringRateAll50ms = np.empty(len(celldb))
laserFiringRateAll50ms = np.empty(len(celldb))
pValLaserAll50ms = np.empty(len(celldb))

for indCell, (indRow, dbRow) in enumerate(celldb.iterrows()):
    oneCell = ephyscore.Cell(dbRow, useModifiedClusters=True)
    
    ephysData, noBehav = oneCell.load('laserPulse')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['laserOn']

    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    # -- For response in first 10 ms --
    baseSpikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 baseRange10ms)
    laserSpikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                  indexLimitsEachTrial,
                                                                  responseRange10ms)
    [testStatistic, pVal] = stats.ranksums(laserSpikeCountMat, baseSpikeCountMat)

    baseFiringRate = np.mean(baseSpikeCountMat)/(baseRange10ms[-1]-baseRange10ms[0])
    laserFiringRate = np.mean(laserSpikeCountMat)/(responseRange10ms[-1]-responseRange10ms[0])

    baseFiringRateAll10ms[indCell] = baseFiringRate
    laserFiringRateAll10ms[indCell] = laserFiringRate
    pValLaserAll10ms[indCell] = pVal
    
    # -- For response in first 50 ms --
    baseSpikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 baseRange50ms)
    laserSpikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                  indexLimitsEachTrial,
                                                                  responseRange50ms)
    [testStatistic, pVal] = stats.ranksums(laserSpikeCountMat, baseSpikeCountMat)

    baseFiringRate = np.mean(baseSpikeCountMat)/(baseRange50ms[-1]-baseRange50ms[0])
    laserFiringRate = np.mean(laserSpikeCountMat)/(responseRange50ms[-1]-responseRange50ms[0])

    baseFiringRateAll50ms[indCell] = baseFiringRate
    laserFiringRateAll50ms[indCell] = laserFiringRate
    pValLaserAll50ms[indCell] = pVal
    
    titleStr = (f'[{indRow}] {oneCell}\tB50={baseFiringRate:0.1f} spk/s ' +
                f'E50={laserFiringRate:0.1f} spk/s (p={pVal:0.3f})')
    print(titleStr)
    
    if 0:
        titleStr = (f'[{indRow}] {oneCell}\nB50={baseFiringRate:0.1f} spk/s ' +
                    f'E50={laserFiringRate:0.1f} spk/s (p={pVal:0.3f})')
        plt.clf()
        plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.')
        plt.title(titleStr)
        plt.show()
        plt.waitforbuttonpress()
    
celldb['baseFiringRate10ms'] = baseFiringRateAll10ms
celldb['laserFiringRate10ms'] = laserFiringRateAll10ms
celldb['pValLaser10ms'] = pValLaserAll10ms
celldb['baseFiringRate50ms'] = baseFiringRateAll50ms
celldb['laserFiringRate50ms'] = laserFiringRateAll50ms
celldb['pValLaser50ms'] = pValLaserAll50ms

celldatabase.save_hdf(celldb, dbPathNew)
