"""
Figure showing learning curves during stage 3.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
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

FIGNAME = 'learning_curve_stage3'
figDataFile = 'fraction_correct_stage3.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)

SAVE_FIGURE = 0
outputDir = '/tmp/'
#avgType = '_no_color' # 'mean' 'median'
avgType = '_labeled' # 'mean' 'median'
figFilename = 'learning_curve_stage3_each_mouse'+avgType # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [9, 7] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeLabels  #figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

np.random.seed(4)

# Exclude dates for very slow mouse and last session stage 3 (failed rigs)
sessionRange = ['2022-02-01', '2022-02-26']
#sessionRange = ['2022-02-01', '2022-02-27'] # Exclude dates for very slow mouse
sessionsToPlot = behavioranalysis.sessions_in_range(sessionRange)

dframe = pd.read_csv(figDataFullPath, index_col=0)
dframe = dframe[sessionsToPlot]

# -- Exclude mice that did not pass criteria 2->3 at the same time --
miceToExclude = ['pamo020'] #[]#['pamo020','pamo024']
dframe = dframe.drop(index=miceToExclude)
#subjects = ['pamo020']
if 0:
    mouseSubset = [f'pamo{s:03d}' for s in range(9,25)]
    mouseSubset = [m for m in mouseSubset if m not in miceToExclude]
    dframe = dframe.loc[mouseSubset]

subjects = list(dframe.index)
nSubjects = len(subjects)

def learningcurve(xval, alpha, beta, lamb, gamma):
    """
    alpha:offset; beta:slope; lamb:max, gamma:min
    """
    return gamma + (lamb-gamma)/(1+np.exp(-(xval-alpha)/beta))

def invlearningcurve(yval, alpha, beta, lamb, gamma):
    """
    Inverse of learningcurve()
    """
    return alpha - beta * np.log((lamb-gamma)/(yval-gamma) - 1)
    
    
plt.clf()
slopes = np.empty(nSubjects)
intercepts = np.empty(nSubjects)
cross = np.empty(nSubjects)
finalperf = np.empty(nSubjects)

for inds,subject in enumerate(subjects):

    perfThisSubject = dframe.loc[subject].to_numpy()
    print(dframe.loc[subject])

    '''
    if 0:
        firstNotNaN = np.flatnonzero(~np.isnan(perfThisSubject))[0]
        perfThisSubject = perfThisSubject[firstNotNaN:]
    else:
        perfThisSubject[np.isnan(perfThisSubject)] = 0.5
    '''
    ########### TEST NaN ############
    #perfThisSubject[19] = np.nan

    
    FIT_TYPE = 'linear'
    #FIT_TYPE = 'sigmoid'
    if FIT_TYPE == 'linear':
        # -- Linear fit --
        dateInds = np.arange(len(perfThisSubject))
        notNaN = ~np.isnan(perfThisSubject)  # Used to mask NaN
        slope, intercept, rval, pval, se = stats.linregress(dateInds[notNaN], perfThisSubject[notNaN])
        #print(f'{subject}: m={slope:0.3f}  b={intercept:0.3f}')
        slopes[inds] = slope
        intercepts[inds] = intercept
        xfit = dateInds[[0,-1]]
        yfit = slope * xfit + intercept
    elif FIT_TYPE == 'sigmoid':
        # -- Sigmoidal fit --
        dateInds = np.arange(len(perfThisSubject))
        p0 = [15.0, 5, 0.8, 0.5]
        try:
            popt, pcov = optimize.curve_fit(learningcurve, dateInds,
                                            perfThisSubject, p0=p0, maxfev=5000)
        except RuntimeError:
            print('Could not fit {subject}')
            popt = p0
        xfit = dateInds
        yfit = learningcurve(xfit, *popt)
        yp0 = learningcurve(xfit, *p0)
        print(popt)
        slopes[inds] = popt[1]
        intercepts[inds] = popt[3]
        cross[inds] = invlearningcurve(70, *popt) 
        finalperf[inds] = learningcurve(24, *popt) # popt[2]
        
    if 1:
        if 0:
            plt.cla()
            plt.plot(dateInds, 100*perfThisSubject,'o', color='0.75', lw=3)
            plt.plot(xfit, 100*yfit, color='b', lw=3)
            plt.axhline(50, ls='--', color='0.8')
            plt.axhline(75, ls='--', color='0.8')
        else:
            plt.cla()
            if subject in studyparams.MICE_ACTIVE_ONLY:
                colorThisSubject = figparams.colors['activeOnly']
            else:
                colorThisSubject = figparams.colors['activePassive']
            plt.axhline(50, ls='--', color='0.8')
            plt.axhline(75, ls='--', color='0.8')
            plt.axhline(70, ls='--', color='0.8')
            plt.plot(dateInds, 100*perfThisSubject,'o-', color='0.85', lw=2, zorder=-1)
            ###plt.plot(dateInds, 100*perfThisSubject,'o-', color=colorThisSubject, lw=2)
            #colorThisSubject = '0.25'
            plt.plot(xfit, 100*yfit, color=colorThisSubject, lw=3)
            #plt.plot(xfit, 100*yp0, color='0.5', lw=3)
           
        plt.yticks(np.arange(50, 110, 10))
        plt.ylim([40, 100])
        plt.ylabel('Correct trials (%)', fontsize=fontSizeLabels)
        plt.xlabel('Sessions in stage 3', fontsize=fontSizeLabels)
        #labelEachCond = [f'{eachCond[indc]} (N={nSubjectsEachCond[indc]})' for indc in [0,1]]
        #plt.legend(lineEachCond[::-1], labelEachCond[::-1], loc='upper left', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        plt.title(f'{subject} - performance stage 3', fontweight='bold')
        plt.show()
        plt.pause(0.01)
        #print(f'Slope ({subject}): {slopes[inds]}')
        print(20*slope + intercept)
        #sys.exit()
        plt.waitforbuttonpress()

#if avgType=='labeled':
#    legend

if FIT_TYPE == 'linear':
    dframe['slope'] = slopes
    dframe['intercept'] = intercepts
    perfThreshold = 0.7
    crossPerfThreshold = (perfThreshold - intercepts)/slopes
    dframe['cross'] = crossPerfThreshold
    finalSessionInd = 25
    dframe['finalPerf'] = slopes*finalSessionInd + intercepts
elif FIT_TYPE == 'sigmoid':
    dframe['slope'] = slopes
    dframe['intercept'] = intercepts
    dframe['cross'] = cross
    dframe['finalPerf'] = finalperf


    
eachCond = ['activeOnly', 'activePassive']
colorEachCond = [figparams.colors['activeOnly'], figparams.colors['activePassive']]

# -- Select animals by performance --
if 1:
    maxSessions = 21
    dframeSubset = dframe[dframe.cross<maxSessions]
    subjects = list(dframeSubset.index)
else:
    dframeSubset = dframe
    
miceActiveOnly = list(set(studyparams.MICE_ACTIVE_ONLY) & set(subjects))
miceActivePassive = list(set(studyparams.MICE_ACTIVE_PASSIVE) & set(subjects))

dataEachCond = [ dframeSubset.loc[miceActiveOnly],
                 dframeSubset.loc[miceActivePassive] ]
nSubjectsEachCond = [len(dataOneCond) for dataOneCond in dataEachCond]

rstat, pval = stats.ranksums(dataEachCond[0].slope, dataEachCond[1].slope)
print(f'Slope: p = {pval:0.4f}')
rstat, pval = stats.ranksums(dataEachCond[0].cross, dataEachCond[1].cross)
print(f'Cross: p = {pval:0.4f}')
rstat, pval = stats.ranksums(dataEachCond[0].finalPerf, dataEachCond[1].finalPerf)
print(f'FinalPerf: p = {pval:0.4f}')

if 0:
    plt.clf()
    plt.subplot(1,3,1)
    for indc, cond in enumerate(eachCond):
        randoffset = 0.1*(2*np.random.rand(nSubjectsEachCond[indc])-1)
        plt.plot(np.tile(indc, nSubjectsEachCond[indc])+randoffset, dataEachCond[indc].slope,
                 'o', mfc='none', mec=colorEachCond[indc], mew=1, ms=6)
    plt.ylabel('Slope')
    plt.xlim([-1,2])
    plt.ylim([0.000, 0.025])
    plt.subplot(1,3,2)
    for indc, cond in enumerate(eachCond):
        randoffset = 0.1*(2*np.random.rand(nSubjectsEachCond[indc])-1)
        plt.plot(np.tile(indc, nSubjectsEachCond[indc])+randoffset, dataEachCond[indc].cross,
                 'o', mfc='none', mec=colorEachCond[indc], mew=1, ms=6)
    plt.ylabel(f'Cross {perfThreshold:0.0%}')
    plt.xlim([-1,2])
    plt.subplot(1,3,3)
    for indc, cond in enumerate(eachCond):
        randoffset = 0.1*(2*np.random.rand(nSubjectsEachCond[indc])-1)
        plt.plot(np.tile(indc, nSubjectsEachCond[indc])+randoffset, dataEachCond[indc].finalPerf,
                 'o', mfc='none', mec=colorEachCond[indc], mew=1, ms=6)
    plt.ylabel(f'Final performance (session {finalSessionInd})')
    plt.xlim([-1,2])
    plt.show()
    
if SAVE_FIGURE:    
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
