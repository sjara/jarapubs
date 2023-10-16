"""
Show distributions of learning speeds for one cohort.
"""

import os
import sys
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


FIGNAME = 'learning_comparison'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

PANELS = [1, 1, 1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_learning_comparison' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 5] # In inches (I'm doubling the size)

fontSizeLabels = figparams.fontSizeLabels + 2
fontSizeTicks = figparams.fontSizeTicks + 2
fontSizePanel = figparams.fontSizePanel

labelPosXtop = [0.01, 0.43, 0.73]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.43, 0.615, 0.815]   # Horiz position for panel labels
labelPosY = [0.95, 0.46]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorEachCond = [figparams.colors['activeOnly'],
                 figparams.colors['activePassive'],
                 figparams.colors['passiveThenActive']]
eachCondLabel = ['A only', 'A + P', 'P : A']


# -- Load distributions --
figDataFile = 'distributions_learning.npz'
figDataFullPath = os.path.join(figDataDir, figDataFile)
distrib = np.load(figDataFullPath, allow_pickle=True)
nComponents = distrib['meanP21'].shape[1] # From the Gaussian mixture model
eachCond = distrib['eachCond']
assert list(eachCond)==['activeOnly', 'activePassive', 'passiveThenActive']
eachCond = ['activeOnly']

#selectedSubjects = studyparams.MICE_ACTIVE_ONLY_COH4
selectedSubjects = studyparams.MICE_ACTIVE_ONLY_ALL
#selectedSubjects = studyparams.MICE_ACTIVE_ONLY_COH2 + studyparams.MICE_ACTIVE_ONLY_COH4
#selectedSubjects = studyparams.MICE_ACTIVE_ONLY_COH3


# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(1, 1)
gsMain.update(left=0.1, right=0.95, top=0.95, bottom=0.2, wspace=0.2, hspace=0.4)
              
# -- Panel: name of panel --
if PANELS[0]:
    #gsMain = gridspec.GridSpec(3, 2)
    #gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.15, wspace=0.3, hspace=0.2)
    #gsDist = gsMain[0,0].subgridspec(3, 2)
    gsDist = gsMain[0,0].subgridspec(1, 2)
    axList = [[0,0,0],[0,0,0]]

    # -- Plot performance at 21 days --
    #histOffset = [ -0.5, -0.5, 0.33 ]
    histOffset = [ -0.6, -0.5, 0.33 ] # -0.6
    for indcond, cond in enumerate(eachCond):
        axList[0][indcond] = plt.subplot(gsDist[indcond, 0])
        indsToPlot = [x in selectedSubjects for x in distrib['subjectsEachCond'][indcond]]
        dataToPlot = distrib['dataP21'][indcond][indsToPlot]
        
        nbins =  16#16
        bins = 100*np.linspace(0.5, 1, nbins) + histOffset[indcond]
        plt.hist(100*dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

        if 0:
            xfit = np.linspace(0.5, 1, 100)
            for indcomp in range(nComponents):
                amp = 4.5*distrib['ampP21'][indcond, indcomp]
                yfit = studyutils.gaussian(xfit, amp, distrib['meanP21'][indcond, indcomp],
                                           distrib['stdP21'][indcond, indcomp], 0)
                plt.plot(100*xfit, yfit, lw=2.5, color=colorEachCond[indcond])
        #plt.ylim([0, 3.5])
        #plt.yticks([0, 2])
        plt.ylabel(f'N mice', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        if indcond==len(eachCond)-1:
            axList[0][indcond].set_xlabel('Performance at 21 days (% correct)', fontsize=fontSizeLabels)
        else:
            axList[0][indcond].set_xticklabels([])
        nMiceLabel = f'N mice = {len(dataToPlot)}'
        inLabelsYpos = 0.8*plt.ylim()[1]
        plt.text(90, inLabelsYpos, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='bold')
        plt.text(90, inLabelsYpos*0.7/0.8, nMiceLabel, fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='normal')

    #axList[0][0].annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
    #                      fontsize=fontSizePanel, fontweight='bold')
    #axList[0][0].annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
    #                      fontsize=fontSizePanel, fontweight='bold')

    # -- Plot time to reach 70% --
    histOffset = [ 0.5, 0.5, 2 ]
    for indcond, cond in enumerate(eachCond):
        axList[0][indcond] = plt.subplot(gsDist[indcond, 1])
        indsToPlot = [x in selectedSubjects for x in distrib['subjectsEachCond'][indcond]]
        dataToPlot = distrib['dataT70'][indcond][indsToPlot]
        nbins = 16
        bins = np.linspace(5, 50, nbins) + histOffset[indcond]
        #bins = np.linspace(5, 150, nbins) + histOffset[indcond]
        plt.hist(dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

        xfit = np.linspace(0, 50, 100)
        for indcomp in range(nComponents):
            amp = 4.5*distrib['ampT70'][indcond, indcomp]
            yfit = studyutils.gaussian(xfit, amp, distrib['meanT70'][indcond, indcomp],
                                       distrib['stdT70'][indcond, indcomp], 0)
            plt.plot(xfit, yfit, lw=2.5, color=colorEachCond[indcond])
        #plt.ylim([0, 3.5])
        #plt.yticks([0, 2])
        #plt.ylabel(f'N mice', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        if indcond==len(eachCond)-1:
            axList[0][indcond].set_xlabel('Time to reach 70% correct (days)', fontsize=fontSizeLabels)
        else:
            axList[0][indcond].set_xticklabels([])
        nMiceLabel = f'N mice = {len(dataToPlot)}'
        inLabelsYpos = 0.8*plt.ylim()[1]
        plt.text(35, inLabelsYpos, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='bold')
        plt.text(35, inLabelsYpos*0.7/0.8, nMiceLabel, fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='normal')
    #axList[0][0].annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
    #                      fontsize=fontSizePanel, fontweight='bold')



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
