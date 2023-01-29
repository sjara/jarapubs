"""
This script provides a template for scripts to generate figures.
Replace this comment with a description of which figure will be plotted by this script.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
#import scipy.stats as stats 
import studyparams
import figparams

FIGNAME = 'figure_name'
figDataFile = 'file_containing_data_for_this_fig.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

PANELS = [1, 1]  # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_figure_name' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.02, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.71]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
laserColor = figparams.colors['blueLaser']

# -- Load data --
figData = np.load(figDataFullPath)

# -- Process data --
pass

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
gsMain = gridspec.GridSpec(1, 2, width_ratios=[0.3, 0.7])
gsMain.update(left=0.1, right=0.98, top=0.95, bottom=0.08, wspace=0.25, hspace=0.3)

# -- Panel: name of panel --
ax0 = plt.subplot(gsMain[0, 0])
ax0.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[0]:
    # Plot stuff
    pass

# -- Panel: name of next panel --
gsPanelsRight = gsMain[0, 1].subgridspec(2, 1, hspace=0.1, height_ratios=[0.2, 0.8])
ax1 = plt.subplot(gsPanelsRight[0, 0])
ax1.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

ax2 = plt.subplot(gsPanelsRight[1, 0])
#ax2.annotate('C', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction',
#             fontsize=fontSizePanel, fontweight='bold')

if PANELS[1]:
    # Plot stuff
    ax1.set_xticklabels([])
    pass

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
