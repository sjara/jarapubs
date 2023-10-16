"""
Find the best number of Gaussians for the Gaussian Mixture Model on performance at 21d.
"""

import sys
sys.path.append('..')
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
reload(studyparams)
from sklearn.mixture import GaussianMixture 

FIGNAME = 'learning_comparison'
figDataFile = 'distributions_learning.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeLabels  #figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

dframe = studyutils.load_stage3(excludeAntibias=1)  # Loads 26 days by default
dframeFit = studyutils.fit_learning_curves(dframe)

subjects = list(dframeFit.index)

eachCond = ['activeOnly', 'activePassive', 'passiveThenActive']
eachCondLabel = ['A only', 'A + P', 'P : A']
colorEachCond = [figparams.colors['activeOnly'], figparams.colors['activePassive'],
                 figparams.colors['passiveThenActive']]
nCond = len(eachCond)

#miceActiveOnly = studyparams.MICE_ACTIVE_ONLY_ALL
#miceActiveOnly = studyparams.MICE_ACTIVE_ONLY_COH4
miceActiveOnly = studyparams.MICE_ACTIVE_ONLY_COH2
miceActivePassive = studyparams.MICE_ACTIVE_PASSIVE
micePassiveThenActive = studyparams.MICE_PASSIVE_THEN_ACTIVE

dataEachCond = [ dframeFit.loc[miceActiveOnly],
                 dframeFit.loc[miceActivePassive],
                 dframeFit.loc[micePassiveThenActive] ]
nSubjectsEachCond = [len(dataOneCond) for dataOneCond in dataEachCond]
subjectsEachCond = [miceActiveOnly, miceActivePassive, micePassiveThenActive]


nCompToTest = [1, 2, 3, 4] # For the Gaussian mixture model
bicP21 = np.empty((len(nCompToTest), nCond))

plt.clf()
gsMain = gridspec.GridSpec(3, len(nCompToTest))
gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.15, wspace=0.3, hspace=0.2)
axList = [[0,0,0],[0,0,0]]

for indnc, nComp in enumerate(nCompToTest):

    if nComp==2:
        # --- For Cohort 2 ---
        meansInitP21 = [ [[0.74], [0.63]], [[0.83], [0.65]] , [[0.85], [0.70]] ]
        weightsInitP21 = [ [0.75, 0.25], [0.67, 0.33] , [0.5, 0.5] ]
    else: #elif nComp==1:
        meansInitP21 = 3*[None]
        weightsInitP21 = 3*[None]

    meanP21 = np.empty((nCond, nComp))
    stdP21 = np.empty((nCond, nComp))
    ampP21 = np.empty((nCond, nComp))
    crossP21 = np.empty(nCond)
    dataP21 = []
    goodLearnersP21 = []

    for indcond, cond in enumerate(eachCond):
        dataToFitP21 = dataEachCond[indcond][f'perfAt21d'].to_numpy()
        dataP21.append(dataToFitP21)
        gmP21 = GaussianMixture(n_components=nComp,
                                weights_init=weightsInitP21[indcond],
                                means_init=meansInitP21[indcond],
                                random_state=None).fit(dataToFitP21[:,np.newaxis])
        meanP21[indcond, :] = gmP21.means_[:,0]
        stdP21[indcond, :] = np.sqrt(gmP21.covariances_[:,0,0])
        ampP21[indcond, :] = gmP21.weights_
        bicP21[indnc, indcond] = gmP21.bic(dataToFitP21[:,np.newaxis])

    if 1:
        histOffset = [ -0.5, -0.5, 0.33 ]
        for indcond, cond in enumerate(eachCond):
            axList[1][indcond] = plt.subplot(gsMain[indcond, indnc])
            dataToPlot = dataEachCond[indcond][f'perfAt21d'].to_numpy()
            #dataToPlot = dataEachCond[indcond][f'perfAt{dayToTestExtra}d'].to_numpy()
            nbins = 16 # 17
            bins = 100*np.linspace(0.5, 1, nbins) + histOffset[indcond]
            plt.hist(100*dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

            xfit = np.linspace(0.5, 1, 100)
            for indcomp in range(nComp):
                amp = 4.5*ampP21[indcond, indcomp]
                yfit = studyutils.gaussian(xfit, amp, meanP21[indcond, indcomp],
                                           stdP21[indcond, indcomp], 0)
                plt.plot(100*xfit, yfit, lw=3, color=colorEachCond[indcond])
            if indcond==0:
                plt.ylim([0, 4.5])
            else:
                plt.ylim([0, 3.5])
            plt.ylabel(f'N mice', fontsize=fontSizeLabels)
            extraplots.boxoff(plt.gca())
            plt.text(50, 2, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                     color=colorEachCond[indcond], fontweight='bold')
        axList[1][0].set_xticklabels([])
        axList[1][1].set_xticklabels([])
        axList[1][2].set_xlabel(f'Performance at 21 days (% correct)', fontsize=fontSizeLabels)
        plt.show()

print(eachCond)        
print(bicP21)

