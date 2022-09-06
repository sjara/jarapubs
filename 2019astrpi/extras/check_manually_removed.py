"""
Double-check that the cells selected by Devin as bad still correspond to bad cells.
"""

import sys
sys.path.append('..')
import os
#import figparams
import studyparams
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import extraplots
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats


#dbPath = '/tmp/testAll.h5'
#dbPath = '/tmp/testAlland1spike.h5'
dbPath = '/tmp/sj_cells_with_soundsessions_20211211.h5'
celldb = celldatabase.load_hdf(dbPath)

toExclude = np.loadtxt(os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME,
                                    'cell_indices_manually_removed.txt'), dtype=int)

gridShape = [4, 6]
nPlots = np.prod(gridShape)
plt.clf()
plt.tight_layout(2.0)
indp = 0
for indc, oneCell in celldb.iterrows():
    axPos = np.unravel_index(indp%nPlots, gridShape)
    plt.subplot2grid(gridShape, axPos)
    plt.cla()
    pColor = 'r' if indc in toExclude else 'k'
    cellStr = (f'{indc} {oneCell.pdepth:0.0f} g{oneCell.egroup}c{oneCell.cluster}')
    plt.yticks([])
    #plt.ylim([-200,200])
    plt.xticks([])
    plt.title(cellStr)
    plt.plot(oneCell.spikeShape, color=pColor, lw=2)
    plt.plot(oneCell.spikeShape-oneCell.spikeShapeSD, color='0.75', lw=1)
    plt.plot(oneCell.spikeShape+oneCell.spikeShapeSD, color='0.75', lw=1)
    if indp%nPlots == (nPlots-1):
        plt.show()
        plt.gcf().suptitle(f'{oneCell.subject} {oneCell.date}', fontweight='bold') 
        input('Press ENTER')
    indp += 1
    




sys.exit()


'''

columns = ['brainArea', 'cluster', 'date', 'egroup', 
           'isiViolations', 'maxDepth', 'nSpikes', 'paradigm', 'pdepth',
           'probe', 'recordingTrack', 'spikePeakAmplitudes',
           'spikePeakTimes', 'spikeShape', 'spikeShapeQuality', 'spikeShapeSD',
           'subject']
celldb = celldatabase.load_hdf(dbPath, columns=columns)
sys.exit()

dbPath = '/tmp/testd1pi048.h5'
cdb = celldatabase.load_hdf(dbPath)
sys.exit()

toExclude = np.loadtxt(os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME,
                                    'cell_indices_manually_removed.txt'), dtype=int)

import h5py
h5file = h5py.File(dbPath, 'r')
for varname, varvalue in h5file['/'].items():
    #if varvalue.dtype.kind == 'S':
    if np.issubdtype(varvalue, np.integer) or np.issubdtype(varvalue, np.floating):
    #if varvalue.dtype == np.object:
        print(varname)

sys.exit()

h5file['ephysTime'][...]
%time x=list(h5file['ephysTime'][...])
%time pd.DataFrame({'x':x})

%time cdb = celldatabase.load_hdf(dbPath, columns=['subject'])


%time x=list(h5file['subject'][...])


behavSuffix 
brainArea
date
ephysTime
info
paradigm
probe
recordingTrack
sessionType
subject
'''
