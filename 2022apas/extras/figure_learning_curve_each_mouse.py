"""
Show example learning curve for one mouse.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pandas as pd
from scipy import stats
from scipy import optimize
from jaratoolbox import settings
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
sys.path.append('..')
import studyparams
import studyutils
import figparams
from importlib import reload
reload(figparams)
reload(studyutils)
reload(studyparams)


SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'extra_all_learning_curves' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [24, 14] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

dframe = studyutils.load_stage3(excludeAntibias=1)  # Loads 26 days by default
#dframe = studyutils.load_stage3(excludeAntibias=0)  # Loads 26 days by default
dframeFit = studyutils.fit_learning_curves(dframe)

if 1:
    groups = [studyparams.MICE_ACTIVE_ONLY_COH2, studyparams.MICE_ACTIVE_PASSIVE,
              studyparams.MICE_PASSIVE_THEN_ACTIVE, studyparams.MICE_ACTIVE_ONLY_COH3]
    colorEachGroup = [figparams.colors['activeOnly'], figparams.colors['activePassive'],
                      figparams.colors['passiveThenActive'], figparams.colors['activeOnly']]
    plt.clf()
    gsMain = gridspec.GridSpec(2, 2)
    #gsMain.update(left=0.05, right=0.98, top=0.95, bottom=0.05, hspace=0.3, wspace=0.3) # Good for 3x3
    gsMain.update(left=0.05, right=0.98, top=0.95, bottom=0.05, hspace=0.2, wspace=0.2)
    subGrid = [3,3]
else:
    groups = [studyparams.MICE_ACTIVE_ONLY_COH4]
    colorEachGroup = [figparams.colors['activeOnly']]
    plt.clf()
    gsMain = gridspec.GridSpec(1, 1)
    subGrid = [4,5]
                      

for indg, group in enumerate(groups):
    dframeGroup = dframe.loc[group]
    gsGroup = gridspec.GridSpecFromSubplotSpec(*subGrid, subplot_spec=gsMain[indg], hspace=0.3, wspace=0.3)
    for inds, (subject, dbrow) in enumerate(dframeGroup.iterrows()):
        perfThisSubject = dbrow.to_numpy()
        dateInds = np.arange(len(perfThisSubject))
        slope, intercept, rval, pval, se = stats.linregress(dateInds, perfThisSubject)
        xfit = dateInds[[0,-1]]
        yfit = dframeFit.slope[subject] * xfit + dframeFit.intercept[subject]

        #thisAx = plt.subplot(*gridShape, inds+1)
        thisAx = plt.subplot(gsGroup[inds])
        plt.axhline(50, ls='--', color='0.9', zorder=-1)
        plt.axhline(70, ls='--', color='0.9', zorder=-1)
        plt.axvline(21, ls='--', color='0.9', zorder=-1)
        plt.plot(dateInds+1, 100*perfThisSubject,'o-', color='0.6', lw=2, zorder=-1)
        plt.plot(xfit+1, 100*yfit, color=colorEachGroup[indg], lw=3)

        '''
        arrowColor = '0.4'
        plt.arrow(15.7, 80, 0, -6, color=arrowColor, width=0.1, head_width=0.5, head_length=2)
        plt.text(17, 82, 'Time to reach 70%', ha='right', color=arrowColor)

        plt.arrow(21, 90-2, 0, -6, color=arrowColor, width=0.1, head_width=0.5, head_length=2)
        plt.text(23, 92-2, 'Performance at 21 days', ha='right', color=arrowColor)
        '''

        plt.yticks(np.arange(50, 110, 10))
        plt.ylim([40, 100])
        plt.ylabel('Correct (%)', fontsize=fontSizeLabels)
        plt.xlabel('Days', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        #plt.title(f'{subject}', fontweight='bold')
        plt.text(2, 93, f'{subject}', color='0', fontweight='bold', ha='left', fontsize=12)
        #plt.show()
        #plt.pause(0.01)
        #sys.exit()
    #plt.tight_layout()
    #plt.waitforbuttonpress()
    plt.pause(0.01)
    
if SAVE_FIGURE:    
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
