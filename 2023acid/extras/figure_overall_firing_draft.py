"""
This script plots comparisons of overall firing rates each reagent.
"""

import os
import sys
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

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_overall_firing' # Do not include extension
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
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldbAll = celldatabase.load_hdf(dbFilename)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
    celldb = celldatabase.load_hdf(dbFilename)
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
    celldb = celldatabase.load_hdf(dbFilename)
else:
    celldb = celldbAll
    
# -- Process data --
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=5)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldbAll, steadyParams, maxChangeFactor)
#posResponse = celldb['preToneGaussianA'] > 0

metrics = ['ToneBaselineFiringRate', 'ToneFiringRateBestFreq', 'ToneAvgEvokedFiringRate']

#           'ToneGaussianA', 'ToneGaussianX0', 'ToneGaussianSigma', 'ToneGaussianY0']

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

for indm, metric in enumerate(metrics):

    thisMetricFiring = np.empty((len(celldb), len(studyparams.REAGENTS)))
    for indr, reagent in enumerate(studyparams.REAGENTS):
        thisMetricFiring[:, indr] = celldb[reagent + metric]

    #thisMetricFiring = thisMetricFiring[responsive & steady & posResponse, :] 
    thisMetricFiring = thisMetricFiring[responsive & steady, :] 
    #thisMetricFiring = thisMetricFiring[responsive, :]  # Control to compare pre-sal to sal-doi

    modIndex1 = studyutils.modulation_index(thisMetricFiring[:,1], thisMetricFiring[:,0])
    modIndex2 = studyutils.modulation_index(thisMetricFiring[:,2], thisMetricFiring[:,1])
    hasModIndex = (~np.isnan(modIndex1) & ~np.isnan(modIndex2))
    wstatMI, pValModInd = stats.wilcoxon(modIndex1[hasModIndex], modIndex2[hasModIndex])
    medianMI1 = np.nanmedian(modIndex1)
    medianMI2 = np.nanmedian(modIndex2)
    print(f'MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')
    
    axs.append(plt.subplot(gsMain[indm]))
    plt.axhline(0, ls='--', color='0.5', lw=1)
    plt.plot([0, 1, 2], thisMetricFiring.T, '-', color='0.75')
    plt.plot([0, 1, 2], thisMetricFiring.T, 'o', mec=pColor, color=pColor)
    for indr, reagent in enumerate(studyparams.REAGENTS):
        medianFR = np.nanmedian(thisMetricFiring[:,indr])
        plt.plot(indr, medianFR, '_', ms=40, color='0.25')
    plt.xlim([-0.5, 2.5])
    plt.ylabel(f'{metric} (spk/s)', fontsize=fontSizeLabels)
    plt.xticks([0, 1, 2], studyparams.REAGENTS, fontsize=fontSizeTicks)
    wstat1, pVal1 = stats.wilcoxon(thisMetricFiring[:,0], thisMetricFiring[:,1])
    wstat2, pVal2 = stats.wilcoxon(thisMetricFiring[:,1], thisMetricFiring[:,2])
    plt.text(0.33, 0.95, f'p = {pVal1:0.3f}', transform=axs[-1].transAxes, ha='center',
             fontsize=fontSizeLabels)
    plt.text(0.66, 0.95, f'p = {pVal2:0.3f}', transform=axs[-1].transAxes, ha='center',
             fontsize=fontSizeLabels)
    #plt.title(f'{cstim} (N={nCellsThisStim})', fontsize=fontSizeLabels)
    #plt.title(f'N = {len(thisMetricFiring)} (p={pValModInd:0.3f})', fontsize=fontSizeLabels)
    plt.title(f'N = {np.sum(hasModIndex)} (p={pValModInd:0.3f})', fontsize=fontSizeLabels)
    plt.show()

sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
