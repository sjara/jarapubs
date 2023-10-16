"""
Plot average learning curves for good and bad mice.
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
figFilename = 'plots_average_learning_curves' # Do not include extension
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
#figDataFile = 'learning_comparison.npz'
figDataFile = 'learning_comparison.20220705.npz'
figDataFullPath = os.path.join(figDataDir, figDataFile)
learningComp = np.load(figDataFullPath, allow_pickle=True)
eachCond = learningComp['eachCond']
assert list(eachCond)==['activeOnly', 'activePassive', 'passiveThenActive']
meanCorrect = learningComp['meanCorrect']
semCorrect = learningComp['semCorrect']
nSubjectsEachCond = learningComp['nSubjectsEachCond']
lastDaysFromDay1 = learningComp['lastDaysFromDay1']
lastDaysStr = f'{lastDaysFromDay1[0]}-{lastDaysFromDay1[-1]}'

'''
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
    print(f'p-value (perf {lastDaysStr} days): {pValPlast[indcomp]:0.04f}')
    print('')
'''

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(1, 2, width_ratios=[0.4, 0.6])
gsMain.update(left=0.07, right=0.99, top=0.95, bottom=0.1, wspace=0.2, hspace=0.4)
              

# -- Panel: learning curves --
ax1 = plt.subplot(gsMain[0, 0])

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


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
