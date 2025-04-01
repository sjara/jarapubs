"""
Estimate SNR of spikes from full voltage traces.
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldb = celldatabase.load_hdf(dbFilename)

spikeShape = np.array(list(celldb.spikeShape))

# -- Process data (copied from figure_overall_firing.py) --
celldbAll = celldb
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=5)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldbAll, steadyParams, maxChangeFactor)
includedCells = responsive & steady

# -- Plot spikes --
fig = plt.gcf()
fig.clf()
#plt.plot(spikeShape[includedCells,:].T)

'''
def plot_spike(spike):
    plt.plot(spike)

extraplots.FlipThrough(plot_spike, spikeShape[includedCells,:], fig=fig)
'''

includedCellsInds = np.where(includedCells)[0]
includedCells = includedCellsInds[[74, 104, 117, 118]]
for ind, spike in enumerate(spikeShape[includedCells,:]):
    plt.cla()
    plt.plot(spike)
    plt.title(f'Cell {ind}')
    plt.ylim(-0.5, 0.5)
    plt.waitforbuttonpress()
plt.show()

# 74, 104, 117, 118  look weird

