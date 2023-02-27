"""
Figure showing the task and example performance for one mouse.
"""

import os
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import studyparams
import studyutils
import figparams
from importlib import reload
reload(figparams)
reload(studyutils)

FIGNAME = 'example_performance'
figDataFile = 'example_performance_{}.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

PANELS = [1, 1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_example_performance' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
#figSize = [15, 6] # In inches (I'm doubling the size)
figSize = [7.5, 2.7] # In inches (I'm doubling the size)

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.38, 0.71]   # Horiz position for panel labels
labelPosY = [0.9, 0.70]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorFit = 'k'
arrowColor = figparams.cp.TangoPalette['Plum1'] #'0.4'

# -- Load data --
subject = 'pamo009'
figData = np.load(figDataFullPath.format(subject))

# -- Processed data --
perfThisSubject = figData['perfEachDayStage3']
dateInds = np.arange(len(perfThisSubject))+1
slope, intercept, rval, pval, se = stats.linregress(dateInds, perfThisSubject)
xfit = dateInds[[0,-1]]
yfit = slope * xfit + intercept
possibleFMslopes = figData['possibleFMslopes']
fractionLeftEachValue = figData['fractionLeftEachValue']
ciLeftEachValue = figData['ciLeftEachValue']
curveParams = figData['curveParams']
             
# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 3, width_ratios=[0.3, 0.35, 0.3], height_ratios=[0.2, 0.8])
gsMain.update(left=0.15, right=0.98, top=0.95, bottom=0.15, wspace=0.4, hspace=0.3)
              
# -- Panel: example learning curve --
ax1 = plt.subplot(gsMain[1, 1])
ax1.annotate('C', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[0]:
    plt.axhline(50, ls='--', color='0.8', lw=0.75)
    plt.axhline(70, ls='--', color='0.8', lw=0.75)
    plt.plot(dateInds, 100*perfThisSubject,'o-', color='0.75', lw=1, ms=4, zorder=-1)
    plt.plot(xfit, 100*yfit, color=colorFit, lw=2)
    plt.arrow(15.7, 80, 0, -6, color=arrowColor, width=0.1, head_width=0.7, head_length=2)
    #plt.text(17, 82, 'Time to reach 70%', ha='right', color=arrowColor)
    plt.text(15.7, 82, 'T to 70%', ha='center', color=arrowColor, fontsize=fontSizeLabels+0)

    plt.arrow(21, 90-2, 0, -6, color=arrowColor, width=0.1, head_width=0.7, head_length=2)
    #plt.text(23, 92-2, 'Performance at 21 days', ha='right', color=arrowColor)
    plt.text(21, 92-2, 'P at 21d', ha='center', color=arrowColor, fontsize=fontSizeLabels+0)

    plt.yticks(np.arange(50, 110, 10))
    plt.ylim([40, 100])
    plt.ylabel('Correct trials (%)', fontsize=fontSizeLabels)
    plt.xlabel('Time after initial shaping (days)', fontsize=fontSizeLabels)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    #plt.text(25.5, 52, f'{subject}', color='0.9', ha='right')


# -- Panel: example psychometric --
ax2 = plt.subplot(gsMain[1, 2])
ax2.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[1]:
    xPad = 0.2 * (possibleFMslopes[-1] - possibleFMslopes[0])
    fitxval = np.linspace(possibleFMslopes[0]-xPad, possibleFMslopes[-1]+xPad, 40)
    fityval = extrastats.psychfun(fitxval, *curveParams)
    xTicks = np.arange(-6, 7, 2)
    hfit = plt.plot(fitxval, 100*fityval, '-', lw=2, color='k')
    (pline, pcaps, pbars, pdots) = studyutils.plot_psychometric(possibleFMslopes,
                                                                fractionLeftEachValue,
                                                                ciLeftEachValue)
    pdots.set_markersize(6)
    pline.set_visible(False)
    plt.xlim([-6, 6])
    plt.ylabel('Leftward choice (%)', fontsize=fontSizeLabels)
    plt.xlabel('FM slope (oct/s)', fontsize=fontSizeLabels)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
