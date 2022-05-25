'''
This script provides a template for scripts to generate figures.
Replace this comment with a description of which figure will be plotted by this script.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
#import scipy.stats as stats 
import figparams

FIGNAME = 'figure_name'
figDataFile = 'file_containing_data_for_this_fig.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

PANELS = [1,1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_figure_name' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.9, 0.48]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
laserColor = figparams.colors['blueLaser']

# -- Load data --
figData = np.load(figDataFullPath)

# -- Processed data --
pass

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 1)
gsMain.update(left=0.15, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)

# -- Panel: name of panel --
ax1 = plt.subplot(gsMain[0, 0])
ax1.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[0]:
    # Plot stuff
    pass

# -- Panel: name of next panel --
ax2 = plt.subplot(gsMain[1, 0])
ax2.annotate('B', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[1]:
    # Plot stuff
    pass

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
