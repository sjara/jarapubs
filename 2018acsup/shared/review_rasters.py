"""
Show example rasters given basic celldb.
"""

import os
import sys
sys.path.append('..')

import pandas as pd
import numpy as np

from jaratoolbox import spikesanalysis
from jaratoolbox import celldatabase
from jaratoolbox import ephyscore

import matplotlib.pyplot as plt


dataDir = '/data/figuresdata/2018acsup/shared/'
dbPath = os.path.join(dataDir,'celldb_2018acsup_photoid_basic.h5')
celldb = celldatabase.load_hdf(dbPath)

for indRow, dbRow in celldb.iterrows():
    oneCell = ephyscore.Cell(dbRow, useModifiedClusters=True)

    #cellInd, dbRow = celldatabase.find_cell(dbase, **cellList[indCell])
    ephysData, noBehav = oneCell.load('laserPulse')

    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['laserOn']
    timeRange = [-0.5, 1]  # In seconds

    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    plt.clf()
    plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.')
    plt.title(str(oneCell))
    plt.show()
    plt.waitforbuttonpress()
    #break


'''
# -- Example cells -- #
cellList = [{'subject' : 'band016',
            'date' : '2016-12-11',
            'depth' : 950,
            'tetrode' : 6,
            'cluster' : 6}, #example AC cell
            ]
'''
