"""
Generate distributions of learning speed.
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
reload(studyparams)
from sklearn.mixture import GaussianMixture 

SAVE_DATA = 0

FIGNAME = 'learning_comparison'
figDataFile = 'distributions_learning.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeLabels  #figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

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


#def gaussian(x, a, x0, sigma, y0):
#    return a*np.exp(-(x-x0)**2/(2*sigma**2))+y0
    

nComp = 2 # For the Gaussian mixture model
if miceActiveOnly == studyparams.MICE_ACTIVE_ONLY_COH2:
    # --- For Cohort 2 ---
    meansInitT70 = [ [[17.7886], [36.2323]], [[12], [25]] , [[12], [25]] ]  #[16.97], [35.25]
    weightsInitT70 = [ [0.73, 0.27], [0.67, 0.33] , [0.67, 0.33] ]
    meansInitP21 = [ [[0.74], [0.63]], [[0.83], [0.65]] , [[0.85], [0.70]] ]
    weightsInitP21 = [ [0.75, 0.25], [0.67, 0.33] , [0.5, 0.5] ]
elif miceActiveOnly == studyparams.MICE_ACTIVE_ONLY_ALL:
    # --- For Cohort 2,3,4 ---
    meansInitT70 = [ [[14], [35]], [[12], [25]] , [[12], [25]] ]  #[16.97], [35.25]
    weightsInitT70 = [ [0.73, 0.27], [0.67, 0.33] , [0.67, 0.33] ]
    meansInitP21 = [ [[0.74], [0.63]], [[0.83], [0.65]] , [[0.85], [0.70]] ]
    weightsInitP21 = [ [0.75, 0.25], [0.67, 0.33] , [0.5, 0.5] ]
# --------------------
meanP21 = np.empty((nCond, nComp))
stdP21 = np.empty((nCond, nComp))
ampP21 = np.empty((nCond, nComp))
crossP21 = np.empty(nCond)
meanT70 = np.empty((nCond, nComp))
stdT70 = np.empty((nCond, nComp))
ampT70 = np.empty((nCond, nComp))
crossT70 = np.empty(nCond)
dataP21 = []
dataT70 = []
goodLearnersP21 = []
goodLearnersT70 = []
for indcond, cond in enumerate(eachCond):
    dataToFitP21 = dataEachCond[indcond][f'perfAt21d'].to_numpy()
    dataP21.append(dataToFitP21)
    gmP21 = GaussianMixture(n_components=nComp,
                            weights_init=weightsInitP21[indcond],
                            means_init=meansInitP21[indcond],
                            random_state=None).fit(dataToFitP21[:,np.newaxis])
    dataToFitT70 = dataEachCond[indcond][f'daysTo70percent'].to_numpy()
    dataT70.append(dataToFitT70)
    gmT70 = GaussianMixture(n_components=nComp,
                            weights_init=weightsInitT70[indcond],
                            means_init=meansInitT70[indcond],
                            random_state=None).fit(dataToFitT70[:,np.newaxis])
    meanP21[indcond, :] = gmP21.means_[:,0]
    stdP21[indcond, :] = np.sqrt(gmP21.covariances_[:,0,0])
    ampP21[indcond, :] = gmP21.weights_
    meanT70[indcond, :] = gmT70.means_[:,0]
    stdT70[indcond, :] = np.sqrt(gmT70.covariances_[:,0,0])
    ampT70[indcond, :] = gmT70.weights_
    # -- Find crossing point between means --
    cross = studyutils.gaussians_mid_cross(meanP21[indcond,0], meanP21[indcond,1],
                                           stdP21[indcond,0], stdP21[indcond,1],
                                           ampP21[indcond, 0],ampP21[indcond, 1] )
    crossP21[indcond] = cross
    cross = studyutils.gaussians_mid_cross(meanT70[indcond,0], meanT70[indcond,1],
                                           stdT70[indcond,0], stdT70[indcond,1],
                                           ampT70[indcond, 0],ampT70[indcond, 1] )
    crossT70[indcond] = cross
    # -- Find good learners --
    goodLearnersP21.append(dataToFitP21 > crossP21[indcond])
    goodLearnersT70.append(dataToFitT70 < crossT70[indcond])
    
    
    
print(f'-- Performance at 21d (means and weights) --')
print(meanP21)
print(ampP21)
print(f'-- Time to 70% (means and weights) --')
print(meanT70)
print(ampT70)

# -- Find good learners --
#subjectsEachCond = [miceActiveOnly, miceActivePassive, micePassiveThenActive]


'''
# -- Find mid-points --
for indcond, cond in enumerate(eachCond):
    xfit = np.linspace(0, 50, 100)
    signFactor = np.sign(meanT70[indcond,1] - meanT70[indcond,0])
    yfitEachCond = []
    for indcomp in range(nComp):
        yfitEachCond.append(studyutils.gaussian(xfit, ampT70[indcond, indcomp],
                                                meanT70[indcond, indcomp],
                                                stdT70[indcond, indcomp], 0))
    diffY = yfitEachCond[1]-yfitEachCond[0]
    plt.clf()
    plt.plot(xfit, diffY)
sys.exit()
'''
    
arrowColor = '0.75'
arrowLabelColor = '0.6'

plt.clf()
gsMain = gridspec.GridSpec(3, 2)
gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.15, wspace=0.3, hspace=0.2)
axList = [[0,0,0],[0,0,0]]

if 1:
    histOffset = [ 0.5, 0.5, 2 ]
    for indcond, cond in enumerate(eachCond):
        axList[0][indcond] = plt.subplot(gsMain[indcond, 1])
        dataToPlot = dataEachCond[indcond][f'daysTo70percent'].to_numpy()
        nbins = 16  # 17
        bins = np.linspace(5, 50, nbins) + histOffset[indcond]
        #bins = np.linspace(-15, 150, 200)
        plt.hist(dataToPlot, bins, color=colorEachCond[indcond], alpha=0.25, rwidth=0.8)

        xfit = np.linspace(0, 50, 100)
        #xfit = np.linspace(10, 150, 100) ## DEBUG
        for indcomp in range(nComp):
            amp = 4.5*ampT70[indcond, indcomp]
            yfit = studyutils.gaussian(xfit, amp, meanT70[indcond, indcomp],
                                       stdT70[indcond, indcomp], 0)
            plt.plot(xfit, yfit, lw=3, color=colorEachCond[indcond])
        if indcond==0:
            plt.ylim([0, 4.5])
        else:
            plt.ylim([0, 3.5])
        plt.ylabel(f'N mice', fontsize=fontSizeLabels)
        plt.text(35, 2, eachCondLabel[indcond], fontsize=1.1*fontSizeLabels, ha='left',
                 color=colorEachCond[indcond], fontweight='bold')
        extraplots.boxoff(plt.gca())
        #sys.exit()
    axList[0][0].set_xticklabels([])
    axList[0][1].set_xticklabels([])
    axList[0][2].set_xlabel(f'Time to reach 70% correct (days)', fontsize=fontSizeLabels)
    #axList[0][0].arrow(34, 3, -16, 0, color=arrowColor, width=0.1, head_width=0.4,
    #          head_length=4, length_includes_head=True)
    #axList[0][0].text(28, 3.3, 'Faster', ha='center', color=arrowLabelColor, fontsize=fontSizeLabels)

    plt.show()

if 1:
    histOffset = [ -0.5, -0.5, 0.33 ]
    for indcond, cond in enumerate(eachCond):
        axList[1][indcond] = plt.subplot(gsMain[indcond, 0])
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
    #plt.arrow(62, 3, 16, 0, color=arrowColor, width=0.1, head_width=0.4,
    #          head_length=4, length_includes_head=True)
    #plt.text(69, 3.3, 'Better', ha='center', color=arrowLabelColor, fontsize=fontSizeLabels)
    plt.show()
    #plt.title(nbins); plt.waitforbuttonpress()


if SAVE_DATA:
    np.savez(figDataFullPath, script=scriptFullPath, nSubjectsEachCond=nSubjectsEachCond,
             subjectsEachCond=subjectsEachCond, eachCond=eachCond,
             meanP21=meanP21, stdP21=stdP21, ampP21=ampP21, crossP21=crossP21,
             meanT70=meanT70, stdT70=stdT70, ampT70=ampT70, crossT70=crossT70,
             dataP21=dataP21, dataT70=dataT70,
             goodLearnersP21=goodLearnersP21, goodLearnersT70=goodLearnersT70)
    print(f'Saved {figDataFullPath}')

