"""
This script plots comparisons of frequency tuning under each reagent.
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
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_freqtuning_R2' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

if len(sys.argv)==2:
    trialSubset = sys.argv[1]
else:
    trialSubset = ''
if trialSubset not in ['', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.02, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.71]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
#colorOddball = figparams.colors['oddball']
#colorStandard = figparams.colors['standard']

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
else:
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Process data --
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldb, steadyParams, maxChangeFactor)
responsive = studyutils.find_tone_responsive_cells(celldb, frThreshold=10)
selective = studyutils.find_freq_selective(celldb, minR2=studyparams.MIN_R_SQUARED)

#metrics = ['ToneBaselineFiringRate', 'ToneFiringRateBestFreq', 'ToneAvgEvokedFiringRate']
metrics = ['ToneGaussianA', 'ToneGaussianX0', 'ToneGaussianSigma', 'ToneGaussianY0']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
#gsMain = gridspec.GridSpec(1, len(studyparams.REAGENTS))
gsMain = gridspec.GridSpec(1, len(metrics))
gsMain.update(left=0.07, right=0.98, top=0.92, bottom=0.08, wspace=0.4, hspace=0.3)
axs = []

pColor = '0.5'

# -- Read frequencies used during tuning --
oneCell = ephyscore.Cell(celldb.iloc[0])
ephysData, bdata = oneCell.load('prePureTones')
stimEachTrial = bdata['currentFreq']
possibleFreq = np.unique(stimEachTrial)
nFreq = len(possibleFreq)
possibleLogFreq = np.log2(possibleFreq)

DEBUG = 0
if DEBUG:
    indRow = 46  # 46 # 55
    indRow = 30
    indRow = 67 # (should fit but doesn't)
    indRow = 435 # inverted, but only inverts on pre
    # 518, 538, 1005, 1006, 1065, 1113, 1125, 1145, 1470 (great fits)
    # 1501 nice fit inverted
    # Check: 1533 raster (very narrow fit
    indRow = 1113 # Good fits to test during running/not
    indRow = 535 # DOI looks flat
    #  [407, 1243, 1285, 1320, 1474, 1480]  These give neg DOI max rate. 1474, 1320 should be positive
    indRow = 1320 # DOI looks flat
    celldbToUse = celldb.iloc[[indRow]]
else:
    celldbToUse = celldb[selective & steady]
    #celldbToUse = celldb

# Cells with decrease: 149, 200, 430, 432, 669, 732, 806, 1145
# Cells with increase: 1113

#celldbToUse = celldb[celldb['preToneGaussianA']<-2000]
#celldbToUse = celldb[celldb['preToneGaussianRsquare'] > studyparams.MIN_R_SQUARED]
#celldbToUse = celldb[(celldb['preToneGaussianSigma']>7) & (celldb['preToneGaussianA']>0)]

for indRow, dbRow in celldbToUse.iterrows():
    #if dbRow['preToneGaussianRsquare'] < studyparams.MIN_R_SQUARED:
    #    pass #continue
    axs = [None, None, None]
    plt.clf()
    for indr, reagent in enumerate(studyparams.REAGENTS):
        fitParams = [dbRow[reagent+'ToneGaussianA'], dbRow[reagent+'ToneGaussianX0'],
                     dbRow[reagent+'ToneGaussianSigma'], dbRow[reagent+'ToneGaussianY0']]
        Rsquared = dbRow[reagent+'ToneGaussianRsquare']
        avgfiringRateEachStim = dbRow[reagent+'ToneFiringRateEachFreq']
        if 0:
            axs[indr] = plt.subplot(1,3,indr+1)
        else:
            axs[indr] = plt.subplot(1,3,indr+1, sharey=axs[0])
        pdots = plt.plot(possibleLogFreq, avgfiringRateEachStim, 'o')
        if not np.isnan(fitParams[0]):
            xvals = np.linspace(possibleLogFreq[0], possibleLogFreq[-1], 60)
            yvals = spikesanalysis.gaussian(xvals, *fitParams)
            pfit = plt.plot(xvals, yvals, '-', lw=3)
        plt.ylabel('Firing rate (Hz)')
        plt.xlabel('Frequency (kHz)')
        xTickLabels = [f'{freq/1000:0.0f}' for freq in possibleFreq]
        plt.xticks(possibleLogFreq, xTickLabels)
        plt.title(f'{reagent} (R2 = {Rsquared:0.4f})\n s={fitParams[2]:0.4f}')
        plt.suptitle(f'{oneCell} [{indRow}]', fontweight='bold')
    plt.draw()
    plt.show()
    plt.waitforbuttonpress()


# -- Strange fits --
# 535
# celldb[selective & posResponse & steady & (celldb['doiToneGaussianSigma']>6)]
