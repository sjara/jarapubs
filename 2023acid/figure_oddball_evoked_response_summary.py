"""
This script plots the evoked responses of each cell under each oddball/standard
condition and under each reagent.
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

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_oddball_summary' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.02, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.71]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorOddball = figparams.colors['oddball']
colorStandard = figparams.colors['standard']

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Process data --
responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=5)


# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
gsMain = gridspec.GridSpec(len(studyparams.REAGENTS), len(studyparams.STIMULI))
gsMain.update(left=0.07, right=0.98, top=0.92, bottom=0.08, wspace=0.4, hspace=0.3)
axs = []

pColor = '0.5'

for indr, reagent in enumerate(studyparams.REAGENTS):
    for inds, stim in enumerate(studyparams.STIMULI):
        axs.append(plt.subplot(gsMain[indr, inds]))
        cstim = stim.capitalize()
        plt.plot([0, 1], [celldb[reagent+cstim+'StandardEvokedFiringRate'][responsive[stim]],
                          celldb[reagent+cstim+'OddballEvokedFiringRate'][responsive[stim]]],
                 '-', color='0.75', lw=1)
        plt.plot(np.zeros(sum(responsive[stim])),
                 celldb[reagent+cstim+'StandardEvokedFiringRate'][responsive[stim]],
                 'o', mec=pColor, color=pColor)
        plt.plot(np.ones(sum(responsive[stim])),
                 celldb[reagent+cstim+'OddballEvokedFiringRate'][responsive[stim]],
                 'o', mfc='w', mew=1.5, color=pColor)
        plt.plot(0, celldb[reagent+cstim+'StandardEvokedFiringRate'][responsive[stim]].median(),
                 '_', ms=20, color='0.25')
        plt.plot(1, celldb[reagent+cstim+'OddballEvokedFiringRate'][responsive[stim]].median(),
                 '_', ms=20, color='0.25')
        plt.xlim([-0.5, 1.5])
        plt.ylabel(f'{reagent}\nEvoked firing rate (spk/s)', fontsize=fontSizeLabels)
        plt.xticks([0, 1], ['Standard', 'Oddball'], fontsize=fontSizeTicks)
        wstat, pVal = stats.wilcoxon(celldb[reagent+cstim+'StandardEvokedFiringRate'][responsive[stim]],
                                     celldb[reagent+cstim+'OddballEvokedFiringRate'][responsive[stim]])

        plt.title(f'{cstim} (N={responsive[stim].sum()})\np = {pVal:0.6f}', fontsize=fontSizeLabels)

    plt.show()
sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
