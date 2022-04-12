"""
Figure showing learning curves during stage 3.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import studyparams
import studyutils
import figparams
from jaratoolbox import settings
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from importlib import reload
reload(figparams)

FIGNAME = 'learning_curve_stage3'
figDataFile = 'fraction_correct_stage3.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

SAVE_FIGURE = 0
outputDir = '/tmp/'
avgType = 'median' # 'mean' 'median'
figFilename = 'learning_curve_stage3_'+avgType # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeLabels  #figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

sessionRange = ['2022-02-01', '2022-02-27']
sessionsToPlot = behavioranalysis.sessions_in_range(sessionRange)

dframe = pd.read_csv(figDataFullPath, index_col=0)
dframe = dframe[sessionsToPlot]

eachCond = ['activeOnly', 'activePassive']
colorEachCond = [figparams.colors['activeOnly'], figparams.colors['activePassive']]

dataEachCond = [ dframe.loc[studyparams.MICE_ACTIVE_ONLY],
                 dframe.loc[studyparams.MICE_ACTIVE_PASSIVE] ]
nSubjectsEachCond = [len(dataOneCond) for dataOneCond in dataEachCond]

plt.clf()
#plt.plot(dframe.median(),'o-', color='k', lw=3); ylim([0.5, 0.9])
plt.axhline(50, ls='--', color='0.8')
lineEachCond = []
if avgType=='mean':
    for indc, cond in enumerate(eachCond):
        avgCorrect = dataEachCond[indc].mean()
        devCorrect = dataEachCond[indc].sem()
        xvals = range(len(avgCorrect))
        plt.fill_between(xvals, 100*(avgCorrect-devCorrect), 100*(avgCorrect+devCorrect),
                         color=colorEachCond[indc], lw=0, alpha=0.25)
        hline, = plt.plot(xvals, 100*avgCorrect, 'o-', color=colorEachCond[indc], lw=3)
        lineEachCond.append(hline)
elif avgType=='median':
    for indc, cond in enumerate(eachCond):
        avgCorrect = dataEachCond[indc].median()
        lowQuantile = dataEachCond[indc].quantile(0.40)
        upQuantile = dataEachCond[indc].quantile(0.60)
        xvals = range(len(avgCorrect))
        #plt.fill_between(xvals, 100*lowQuantile, 100*upQuantile,
        #                 color=colorEachCond[indc], alpha=0.25)
        hline, = plt.plot(xvals, 100*avgCorrect, 'o-', color=colorEachCond[indc], lw=3)
        lineEachCond.append(hline)
    #plt.plot(100*dfAO.median(),'o-', color=colorAO, lw=3);
    #plt.plot(100*dfAP.median(),'o-', color=colorAP, lw=3);
#plt.xticks(rotation=80)
plt.yticks(np.arange(50, 110, 10))
plt.ylim([43, 90])
plt.ylabel('Correct trials (%)', fontsize=fontSizeLabels)
plt.xlabel('Sessions in stage 3', fontsize=fontSizeLabels)
labelEachCond = [f'{eachCond[indc]} (N={nSubjectsEachCond[indc]})' for indc in [0,1]]
plt.legend(lineEachCond[::-1], labelEachCond[::-1], loc='upper left', fontsize=fontSizeLabels)
extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
extraplots.boxoff(plt.gca())
plt.title(f'Average performance across mice ({avgType})', fontweight='bold')
plt.show()

if SAVE_FIGURE:    
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
