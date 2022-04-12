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

FIGNAME = 'learning_curve_stage4'
figDataFile = 'measurements_stage4.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'learning_curve_stage4' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

###sessionRange = ['2022-02-28', '2022-03-20']
###sessionsToPlot = behavioranalysis.sessions_in_range(sessionRange)

dfCurveParams = pd.read_csv(figDataFullPath, index_col=[0,1])

#dfCurveParams = dfCurveParams.drop('pamo019', level=1); print('*** Dropping pamo019 ***')

###dframe = dframe[sessionsToPlot]
#dframe = 1/dfCurveParams.xs('invSlope')
#dframe = dfCurveParams.xs('bias')
#dframe = 1-dfCurveParams.xs('lapseHigh')-dfCurveParams.xs('lapseLow')
dframe = dfCurveParams.xs('nCorrect')/dfCurveParams.xs('nValid')

eachCond = ['activeOnly', 'activePassive']
colorEachCond = [figparams.colors['activeOnly'], figparams.colors['activePassive']]

dataEachCond = [ dframe.loc[studyparams.MICE_ACTIVE_ONLY],
                 dframe.loc[studyparams.MICE_ACTIVE_PASSIVE] ]

nSubjectsEachCond = [len(dataOneCond) for dataOneCond in dataEachCond]

plt.clf()
#plt.plot(dframe.median(),'o-', color='k', lw=3); ylim([0.5, 0.9])
#plt.axhline(50, ls='--', color='0.8')
lineEachCond = []
if 1:
    for indc, cond in enumerate(eachCond):
        avgParam = dataEachCond[indc].mean()
        devParam = dataEachCond[indc].sem()
        xvals = range(len(avgParam))
        plt.fill_between(xvals, avgParam-devParam, avgParam+devParam,
                         color=colorEachCond[indc], lw=0, alpha=0.25)
        hline, = plt.plot(xvals, avgParam, 'o-', color=colorEachCond[indc], lw=3)
        lineEachCond.append(hline)
else:
    for indc, cond in enumerate(eachCond):
        avgParam = dataEachCond[indc].median()
        lowQuantile = dataEachCond[indc].quantile(0.40)
        upQuantile = dataEachCond[indc].quantile(0.60)
        xvals = range(len(avgParam))
        #plt.fill_between(xvals, 100*lowQuantile, 100*upQuantile,
        #                 color=colorEachCond[indc], alpha=0.25)
        hline, = plt.plot(xvals, avgParam, 'o-', color=colorEachCond[indc], lw=3)
        lineEachCond.append(hline)
    #plt.plot(100*dfAO.median(),'o-', color=colorAO, lw=3);
    #plt.plot(100*dfAP.median(),'o-', color=colorAP, lw=3);
#plt.xticks(rotation=80)
#plt.yticks(np.arange(50, 110, 10))
#plt.ylim([43, 90])
plt.ylabel('Some param', fontsize=fontSizeLabels)
plt.xlabel('Sessions in stage 4', fontsize=fontSizeLabels)
labelEachCond = [f'{eachCond[indc]} (N={nSubjectsEachCond[indc]})' for indc in [0,1]]
plt.legend(lineEachCond[::-1], labelEachCond[::-1], loc='upper left', fontsize=fontSizeLabels)
extraplots.boxoff(plt.gca())
plt.show()

if SAVE_FIGURE:    
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
