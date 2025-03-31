"""
Compare responses during running and not-running trials.
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


# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
dbFilenameRun = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
dbFilenameNotrun = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
condLabels = ['notrunning', 'running']
celldbAll = celldatabase.load_hdf(dbFilename)
celldbs = [celldatabase.load_hdf(dbFilenameNotrun),
           celldatabase.load_hdf(dbFilenameRun)]

# -- Process data --
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=10, allreagents=True)
steady = studyutils.find_tone_steady_cells(celldbAll, maxChangeFactor=1.2)

metrics = ['ToneBaselineFiringRate', 'ToneFiringRateBestFreq', 'ToneAvgEvokedFiringRate']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
gsMain = gridspec.GridSpec(len(studyparams.REAGENTS), len(metrics))
gsMain.update(left=0.07, right=0.98, top=0.92, bottom=0.08, wspace=0.4, hspace=0.3)
axs = []

pColor = '0.5'



for indr, reagent in enumerate(studyparams.REAGENTS):
    #baselineFiring = np.empty((len(celldbAll), len(condLabels)))
    #evokedFiring = np.empty((len(celldbAll), len(condLabels)))

    for indm, metric in enumerate(metrics):
        thisMetricFiring = np.empty((len(celldbAll), len(condLabels)))
        for indc, condLabel in enumerate(condLabels):
            thisMetricFiring[:, indc] = celldbs[indc][reagent + metric]
        #thisMetricFiring = thisMetricFiring[responsive & steady, :] 
        thisMetricFiring = thisMetricFiring[responsive, :]  # Control

        axs.append(plt.subplot(gsMain[indr, indm]))
        plt.axhline(0, ls='--', color='0.5', lw=1)
        plt.plot([0, 1], thisMetricFiring.T, '-', color='0.75')
        plt.plot([0, 1], thisMetricFiring.T, 'o', mec=pColor, color=pColor)
        for indc, condLabel in enumerate(condLabels):
            medianFR = np.nanmedian(thisMetricFiring[:,indc])
            print(f'{reagent} {metric} {condLabel} median: {medianFR}')
            plt.plot(indc, medianFR, '_', ms=40, color='0.25')
        plt.xlim([-0.5, 1.5])
        plt.ylabel(f'{metric} (spk/s)', fontsize=fontSizeLabels)
        plt.xticks([0, 1], condLabels, fontsize=fontSizeTicks)
        wstat1, pVal1 = stats.wilcoxon(thisMetricFiring[:,0], thisMetricFiring[:,1])
        #wstat2, pVal2 = stats.wilcoxon(thisMetricFiring[:,1], thisMetricFiring[:,2])
        plt.text(0.5, 0.95, f'p = {pVal1:0.3f}', transform=axs[-1].transAxes, ha='center',
                 fontsize=fontSizeLabels)
        #plt.text(0.66, 0.95, f'p = {pVal2:0.3f}', transform=axs[-1].transAxes, ha='center',
        #         fontsize=fontSizeLabels)
        #plt.title(f'{cstim} (N={nCellsThisStim})', fontsize=fontSizeLabels)
        plt.title(f'{reagent}  (N = {len(thisMetricFiring)})', fontsize=fontSizeLabels)
        plt.show()

sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
