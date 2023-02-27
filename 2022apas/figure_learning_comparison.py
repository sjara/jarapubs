"""
Figure comparing learning across cohorts.
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


learnerGroup = 'fast'  # 'slow'
#learnerGroup = 'slow'
FIGNAME = 'learning_comparison'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

PANELS = [1, 1, 1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_learning_comparison' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
#figSize = [15, 6] # In inches (I'm doubling the size)
figSize = [7.5, 5] # In inches (I'm doubling the size)

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosXtop = [0.01, 0.43, 0.73]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.43, 0.615, 0.815]   # Horiz position for panel labels
labelPosY = [0.95, 0.46]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorEachCond = [figparams.colors['activeOnly'],
                 figparams.colors['activePassive'],
                 figparams.colors['passiveThenActive']]
eachCondLabel = ['A only', 'A + P', 'P : A']

# -- Load learning curves --
#figDataFile = 'average_learning_curves.npz'
figDataFile = f'learning_comparison_{learnerGroup}.npz'
figDataFullPath = os.path.join(figDataDir, figDataFile)
learningComp = np.load(figDataFullPath, allow_pickle=True)
eachCond = learningComp['eachCond']
assert list(eachCond)==['activeOnly', 'activePassive', 'passiveThenActive']
meanCorrect = learningComp['meanCorrect']
semCorrect = learningComp['semCorrect']
nSubjectsEachCond = learningComp['nSubjectsEachCond']
lastDaysFromDay1 = learningComp['lastDaysFromDay1']
lastDaysStr = f'{lastDaysFromDay1[0]}-{lastDaysFromDay1[-1]}'

# -- Load distributions --
figDataFile = 'distributions_learning.npz'
figDataFullPath = os.path.join(figDataDir, figDataFile)
distrib = np.load(figDataFullPath, allow_pickle=True)
nComponents = distrib['meanP21'].shape[1] # From the Gaussian mixture model

# -- Process data --
possibleComp = [[0,1], [0,2], [1,2]] # Possible comparisons
nComparisons = len(possibleComp)
pValT70 = np.empty(nComparisons)
pValP21 = np.empty(nComparisons)
pValPlast = np.empty(nComparisons)
for indcomp, thisComp in enumerate(possibleComp):
    print(f'Comparing: {eachCond[thisComp[0]]} to {eachCond[thisComp[1]]}')
    wstat, pValT70[indcomp] = stats.ranksums(learningComp['dataT70'][thisComp[0]],
                                             learningComp['dataT70'][thisComp[1]])
    print(f'p-value (time to 70%): {pValT70[indcomp]:0.04f}')
    wstat, pValP21[indcomp] = stats.ranksums(learningComp['dataP21'][thisComp[0]],
                                             learningComp['dataP21'][thisComp[1]])
    print(f'p-value (perf at 21d): {pValP21[indcomp]:0.04f}')
    wstat, pValPlast[indcomp] = stats.ranksums(learningComp['dataPlast'][thisComp[0]],
                                               learningComp['dataPlast'][thisComp[1]])
    meanPlast0 = np.mean(learningComp['dataPlast'][thisComp[0]])
    meanPlast1 = np.mean(learningComp['dataPlast'][thisComp[1]])
    print(f'p-value (perf {lastDaysStr} days): {pValPlast[indcomp]:0.04f}' +
          f"  ({meanPlast0:0.1%} vs {meanPlast1:0.1%})")
    print('')


# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 2, width_ratios=[0.4, 0.6], height_ratios=[0.5, 0.5])
gsMain.update(left=0.07, right=0.99, top=0.95, bottom=0.1, wspace=0.2, hspace=0.4)
              
# -- Panel: name of panel --
if PANELS[0]:
    #gsMain = gridspec.GridSpec(3, 2)
    #gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.15, wspace=0.3, hspace=0.2)
    gsDist = gsMain[0,1].subgridspec(3, 2)
    axList = [[0,0,0],[0,0,0]]

    # -- Plot performance at 21 days --
    histOffset = [ -0.5, -0.5, 0.33 ]
    for indcond, cond in enumerate(eachCond):
        axList[0][indcond] = plt.subplot(gsDist[indcond, 0])
        dataToPlot = distrib['dataP21'][indcond]
        nbins = 16
        bins = 100*np.linspace(0.5, 1, nbins) + histOffset[indcond]
        plt.hist(100*dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

        xfit = np.linspace(0.5, 1, 100)
        for indcomp in range(nComponents):
            amp = 4.5*distrib['ampP21'][indcond, indcomp]
            yfit = studyutils.gaussian(xfit, amp, distrib['meanP21'][indcond, indcomp],
                                       distrib['stdP21'][indcond, indcomp], 0)
            plt.plot(100*xfit, yfit, lw=2.5, color=colorEachCond[indcond])
        plt.ylim([0, 3.5])
        plt.yticks([0, 2])
        plt.ylabel(f'N mice', fontsize=fontSizeLabels)
        plt.text(50, 2, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='bold')
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
    axList[0][0].set_xticklabels([])
    axList[0][1].set_xticklabels([])
    axList[0][2].set_xlabel('Performance at 21 days (% correct)', fontsize=fontSizeLabels)

    axList[0][0].annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
                          fontsize=fontSizePanel, fontweight='bold')
    axList[0][0].annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
                          fontsize=fontSizePanel, fontweight='bold')

    # -- Plot time to reach 70% --
    histOffset = [ 0.5, 0.5, 2 ]
    for indcond, cond in enumerate(eachCond):
        axList[0][indcond] = plt.subplot(gsDist[indcond, 1])
        dataToPlot = distrib['dataT70'][indcond]
        nbins = 16
        bins = np.linspace(5, 50, nbins) + histOffset[indcond]
        plt.hist(dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

        xfit = np.linspace(0, 50, 100)
        for indcomp in range(nComponents):
            amp = 4.5*distrib['ampT70'][indcond, indcomp]
            yfit = studyutils.gaussian(xfit, amp, distrib['meanT70'][indcond, indcomp],
                                       distrib['stdT70'][indcond, indcomp], 0)
            plt.plot(xfit, yfit, lw=2.5, color=colorEachCond[indcond])
        plt.ylim([0, 3.5])
        plt.yticks([0, 2])
        #plt.ylabel(f'N mice', fontsize=fontSizeLabels)
        plt.text(35, 2, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='bold')
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
    axList[0][0].set_xticklabels([])
    axList[0][1].set_xticklabels([])
    axList[0][2].set_xlabel('Time to reach 70% correct (days)', fontsize=fontSizeLabels)
    
    axList[0][0].annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
                          fontsize=fontSizePanel, fontweight='bold')


# -- Panel: learning curves --
ax1 = plt.subplot(gsMain[1, 0])
ax1.annotate('D', xy=(labelPosXbot[0],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

if PANELS[1]:
    plt.axhline(50, ls='--', color='0.8')
    lineEachCond = []
    for indc, cond in enumerate(eachCond):
        xvals = range(1, len(meanCorrect[indc])+1)
        plt.fill_between(xvals, 100*(meanCorrect[indc]-semCorrect[indc]),
                         100*(meanCorrect[indc]+semCorrect[indc]),
                         color=colorEachCond[indc], lw=0, alpha=0.25)
        hline, = plt.plot(xvals, 100*meanCorrect[indc], 'o-', color=colorEachCond[indc],
                          lw=2, ms=4)
        lineEachCond.append(hline)
    plt.yticks(np.arange(50, 110, 10))
    plt.ylim([40, 100])
    plt.ylabel('Correct trials (%)', fontsize=fontSizeLabels)
    plt.xlabel('Time after initial shaping (days)', fontsize=fontSizeLabels)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    legendLabels = [f'{lab} ({n} mice)' for lab,n in zip(eachCondLabel,nSubjectsEachCond)]
    plt.legend(lineEachCond, legendLabels, loc='upper left', fontsize=fontSizeLabels, frameon=False)

# -- Panel: comparisons across conditions --
if PANELS[2]:
    markerSizeComp = 4
    medianWidth = 0.3
    gsComp = gsMain[1,1].subgridspec(1, 3, wspace=0.9)
    axCompList = [0, 0, 0]

    # -- Plot time to reach 70% --
    axCompList[0] = plt.subplot(gsComp[0, 0])
    for indcond, cond in enumerate(eachCond):
        dataToPlot = learningComp['dataT70'][indcond]
        medianThisCond = np.median(dataToPlot)
        offset = extraplots.spread_offsets(dataToPlot, 0.1, 0.2)
        plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                 'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
        plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=2,
                 color=colorEachCond[indcond])
    plt.ylabel(f'Time to reach 70% correct (days)', fontsize=fontSizeLabels)
    plt.ylim([0, 30])
    plt.xlim([-0.75, 2.75])
    plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
    if pValT70[0]<0.05:
        hs, hl = extraplots.significance_stars([0, 1], 27, yLength=0.7, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValT70[1]<0.05:
        hs, hl = extraplots.significance_stars([0, 2], 29, yLength=0.7, gapFactor=0.2, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValT70[2]<0.05:
        hs, hl = extraplots.significance_stars([1, 2], 26, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    #if pVal70pc1<0.05:
    #    extraplots.significance_stars([0, 1], 26, yLength=0.5, gapFactor=0.2)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    axCompList[0].annotate('E', xy=(labelPosXbot[1],labelPosY[1]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')

    # -- Performance at 21 days --
    axCompList[1] = plt.subplot(gsComp[0, 1])
    for indcond, cond in enumerate(eachCond):
        dataToPlot = learningComp['dataP21'][indcond]
        medianThisCond = np.median(dataToPlot)
        offset = extraplots.spread_offsets(dataToPlot, 0.1, 0.005)
        plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, 100*dataToPlot,
                 'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
        plt.plot(indcond+medianWidth*np.array([-1,1]), 100*np.tile(medianThisCond,2), lw=2,
                 color=colorEachCond[indcond])
    plt.ylabel(f'Estimated performance\nat 21 days (% correct)', fontsize=fontSizeLabels)
    plt.ylim([50, 100])
    plt.xlim([-0.75, 2.75])
    plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
    if pValP21[0]<0.05:
        hs, hl = extraplots.significance_stars([0, 1], 96, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValP21[1]<0.05:
        hs, hl = extraplots.significance_stars([0, 2], 99, yLength=0.9, gapFactor=0.2, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValP21[2]<0.05:
        hs, hl = extraplots.significance_stars([1, 2], 96, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    axCompList[1].annotate('F', xy=(labelPosXbot[2],labelPosY[1]), xycoords='figure fraction',
                           fontsize=fontSizePanel, fontweight='bold')

    # -- Actual performance last 4 days --
    axCompList[2] = plt.subplot(gsComp[0, 2])
    for indcond, cond in enumerate(eachCond):
        dataToPlot = learningComp['dataPlast'][indcond]
        medianThisCond = np.median(dataToPlot)
        offset = extraplots.spread_offsets(dataToPlot, 0.1, 0.005)
        plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, 100*dataToPlot,
                 'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
        plt.plot(indcond+medianWidth*np.array([-1,1]), 100*np.tile(medianThisCond,2), lw=2,
                 color=colorEachCond[indcond])
    plt.ylabel(f'Actual performance\nat {lastDaysStr} days (% correct)', fontsize=fontSizeLabels)
    plt.ylim([50, 100])
    plt.xlim([-0.75, 2.75])
    plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
    if pValPlast[0]<0.05:
        hs, hl = extraplots.significance_stars([0, 1], 96, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValPlast[1]<0.05:
        hs, hl = extraplots.significance_stars([0, 2], 99, yLength=0.9, gapFactor=0.2, starSize=6)
        plt.setp(hl, lw=0.75)
    if pValPlast[2]<0.05:
        hs, hl = extraplots.significance_stars([1, 2], 96, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    axCompList[2].annotate('G', xy=(labelPosXbot[3],labelPosY[1]), xycoords='figure fraction',
                           fontsize=fontSizePanel, fontweight='bold')

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
