"""
Generate data for comparison of learning speed across conditions.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy import stats
from scipy import optimize
import studyparams
import studyutils
import figparams
from jaratoolbox import settings
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from importlib import reload
reload(figparams)
reload(studyutils)

SAVE_DATA = 0

learnerGroup = 'fast'
#learnerGroup = 'slow'
FIGNAME = 'learning_comparison'
figDataFile = f'learning_comparison_{learnerGroup}.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

# -- Optional --
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()

dframe = studyutils.load_stage3(excludeAntibias=1)  # Loads 26 days by default
dframeFit = studyutils.fit_learning_curves(dframe)

subjects = list(dframeFit.index)

eachCond = ['activeOnly', 'activePassive', 'passiveThenActive']
eachCondLabel = ['A only', 'A + P', 'P : A']
colorEachCond = [figparams.colors['activeOnly'], figparams.colors['activePassive'],
                 figparams.colors['passiveThenActive']]
nCond = len(eachCond)

'''
fastSubjects, slowSubjects = studyutils.find_fast_learners(dframeFit=dframeFit)
#subjectsToInclude = fastSubjects + slowSubjects
subjectsToInclude = fastSubjects
#subjectsToInclude = slowSubjects
'''
fastSubjectsEachCond, eachCond = studyutils.mice_each_condition(learnerGroup, aoc=[2])
#fastSubjectsEachCond, eachCond = studyutils.mice_each_condition('fast', aoc=[2])
#fastSubjectsEachCond, eachCond = studyutils.mice_each_condition('slow', aoc=[2])
#fastSubjectsEachCond, eachCond = studyutils.mice_each_condition('fast', aoc=[2,3,4])
#fastSubjectsEachCond, eachCond = studyutils.mice_each_condition('all')   ### DEBUG
assert eachCond==['activeOnly', 'activePassive', 'passiveThenActive']
miceEachCond = fastSubjectsEachCond

'''
#miceActiveOnly = list(set(studyparams.MICE_ACTIVE_ONLY_ALL) & set(subjectsToInclude))
miceActiveOnly = list(set(studyparams.MICE_ACTIVE_ONLY_COH2) & set(subjectsToInclude))
miceActivePassive = list(set(studyparams.MICE_ACTIVE_PASSIVE) & set(subjectsToInclude))
micePassiveThenActive = list(set(studyparams.MICE_PASSIVE_THEN_ACTIVE) & set(subjectsToInclude))

print('WARNING! excluding one P:A mouse, until I figure out the boundary for each condition')
micePassiveThenActive.remove('pamo035')

# -- Calculate average learning curves --
dataEachCond = [ dframe.loc[miceActiveOnly],
                 dframe.loc[miceActivePassive],
                 dframe.loc[micePassiveThenActive] ]
'''

# -- Calculate average learning curves and store data for each mouse --
dataEachCond = [dframe.loc[miceEachCond[0],:],
                dframe.loc[miceEachCond[1],:],
                dframe.loc[miceEachCond[2],:]]
dataFitEachCond = [ dframeFit.loc[miceEachCond[0]],
                    dframeFit.loc[miceEachCond[1]],
                    dframeFit.loc[miceEachCond[2]] ]

nCond = len(eachCond)
nSubjectsEachCond = [len(dataOneCond) for dataOneCond in dataEachCond]
nSessionEachCond = [dataOneCond.shape[1] for dataOneCond in dataEachCond]
assert np.all(np.array(nSessionEachCond)==nSessionEachCond[0]) # Check same number of sessions
nSessions = nSessionEachCond[0]

avgType = 'mean'
#avgType = 'median' # For testing

lastDaysFromDay1 = [23, 24, 25, 26] # Days to estimate performance at the end of S3 (days start at 1)
lastDaysRange = np.array(lastDaysFromDay1)-1
dframeAvgLastSessions = dframe.iloc[:,lastDaysRange].mean(axis=1)

meanCorrect = np.empty((nCond,nSessions))
semCorrect = np.empty((nCond,nSessions))
dataP21 = []
dataT70 = []
dataPlast = []
for indcond, cond in enumerate(eachCond):
    if avgType=='mean':
        meanCorrect[indcond,:] = dataEachCond[indcond].mean()
        semCorrect[indcond,:] = dataEachCond[indcond].sem()
    elif avgType=='median':
        meanCorrect[indcond,:] = dataEachCond[indcond].median()
        semCorrect[indcond,:] = np.zeros(meanCorrect[indcond,:].shape) # TEMPORARY
    dataToFitP21 = dataFitEachCond[indcond][f'perfAt21d'].to_numpy()
    dataP21.append(dataToFitP21)
    dataToFitT70 = dataFitEachCond[indcond][f'daysTo70percent'].to_numpy()
    dataT70.append(dataToFitT70)
    dataThisCondPlast = dframeAvgLastSessions.loc[miceEachCond[indcond]]
    dataPlast.append(dataThisCondPlast)

'''        
# -- Calculate measurements for comparisons across conditions --
dataFitEachCond = [ dframeFit.loc[miceActiveOnly],
                    dframeFit.loc[miceActivePassive],
                    dframeFit.loc[micePassiveThenActive] ]
lateRange = np.arange(22,26)
dframeAvgLastSessions = dframe.iloc[:,lateRange].mean(axis=1)
lastPerfEachCond = [ dframeAvgLastSessions.loc[miceActiveOnly],
                     dframeAvgLastSessions.loc[miceActivePassive],
                     dframeAvgLastSessions.loc[micePassiveThenActive] ]
'''

if 1:
    plt.clf()
    gsMain = gridspec.GridSpec(1, 1)
    gsMain.update(left=0.13, right=0.98, top=0.9, bottom=0.15)
    ax0 = plt.subplot(gsMain[0, 0])

    fontSizeLabels = figparams.fontSizeLabels
    fontSizeTicks = figparams.fontSizeLabels  #figparams.fontSizeTicks
    fontSizePanel = figparams.fontSizePanel
    
    plt.axhline(50, ls='--', color='0.8')
    lineEachCond = []
    for indc, cond in enumerate(eachCond):
        xvals = range(1, len(meanCorrect[indc])+1)
        plt.fill_between(xvals, 100*(meanCorrect[indc]-semCorrect[indc]),
                         100*(meanCorrect[indc]+semCorrect[indc]),
                         color=colorEachCond[indc], lw=0, alpha=0.25)
        hline, = plt.plot(xvals, 100*meanCorrect[indc], 'o-', color=colorEachCond[indc], lw=3)
        lineEachCond.append(hline)
    plt.yticks(np.arange(50, 110, 10))
    plt.ylim([40, 100])
    plt.ylabel('Correct trials (%)', fontsize=fontSizeLabels)
    plt.xlabel('Time after initial shaping (days)', fontsize=fontSizeLabels)
    extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
    extraplots.boxoff(plt.gca())
    legendLabels = [f'{lab} (N = {n} mice)' for lab,n in zip(eachCondLabel,nSubjectsEachCond)]
    plt.legend(lineEachCond[::1], legendLabels[::1], loc='upper left', fontsize=fontSizeLabels)
    plt.title(f'Learning performance, average across mice', fontweight='bold')
    plt.show()

if SAVE_DATA:
    np.savez(figDataFullPath, script=scriptFullPath, nSubjectsEachCond=nSubjectsEachCond,
             eachCond=eachCond, meanCorrect=meanCorrect, semCorrect=semCorrect,
             dataP21=dataP21, dataT70=dataT70, dataPlast=dataPlast,
             lastDaysFromDay1=lastDaysFromDay1)
    print(f'Saved {figDataFullPath}')
