"""
This script plots comparisons of trial-to-trial variabilty
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_variability' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [14, 6.5] # In inches

if len(sys.argv)==2:
    trialSubset = sys.argv[1]
else:
    trialSubset = ''
if trialSubset not in ['', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.01, 0.27, 0.52, 0.77]   # Horiz position for panel labels
labelPosY = [0.95, 0.48]    # Vert position for panel labels

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_variability.h5')
celldbAll = celldatabase.load_hdf(dbFilename)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_variability_running.h5')
celldbRunning = celldatabase.load_hdf(dbFilename)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_variability_notrunning.h5')
celldbNotRunning = celldatabase.load_hdf(dbFilename)
celldbs = [celldbAll, celldbNotRunning, celldbRunning]
runningConditions = ['all', 'notRunning', 'running']
runningCondLabels = ['All trials', 'Not running', 'Running']

# -- Calculate Fano factor --
roundPossibleFreq = np.round(np.logspace(np.log10(2000), np.log10(40000), 16)) # HARDCODED
for celldb in celldbs:
    roundBestFreq = np.round(celldb['salineToneBestFreq'])
    for indr, reagent in enumerate(studyparams.REAGENTS):
        fanoFactorEachStim = np.array(list(celldb[reagent+'ToneTrialToTrialVarSpks'] /
                                           celldb[reagent+'ToneTrialToTrialMeanSpks']))
        celldb[reagent+'FanoFactorAvg'] = np.nanmean(fanoFactorEachStim, axis=1)
        bestFreqInd = [np.where(roundPossibleFreq==freq)[0][0] for freq in roundBestFreq]
        fanoFactorBestFreq = fanoFactorEachStim[np.arange(len(fanoFactorEachStim)), bestFreqInd]
        celldb[reagent+'FanoFactorBestFreq'] = np.array(fanoFactorBestFreq)

# -- Process data --
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=5)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldbAll, steadyParams, maxChangeFactor)
#posResponse = celldb['preToneGaussianA'] > 0
#metrics = ['ToneBaselineFiringRate', 'ToneFiringRateBestFreq', 'ToneAvgEvokedFiringRate']
#metrics = ['FanoFactorAvg', 'ToneAvgEvokedFiringRate']
#metricLabels = ['Fano factor', 'Avg evoked firing']
metrics = ['FanoFactorAvg', 'FanoFactorBestFreq']
metricLabels = ['Fano factor Avg', 'Fano factor BF']
metricLims = [[0, 2], [0, 2]]
metricTicks = [[0, 20, 40, 60], [0, 40, 80]]

# -- Estimate how many trials in each session under each running condition. --
nTrialsEachRunningCond = [[], []]
for inds, subject in enumerate(studyparams.SUBJECTS):
    for indc, celldbThisCond in enumerate(celldbs[1:]):
        celldbThisSubject = celldbThisCond.query('subject==@subject')
        datesThisSubject = list(celldbThisSubject.date.unique())
        nTrialsThisRunning = []
        for oneDate in datesThisSubject:
            celldbThisDate = celldbThisSubject.query('date==@oneDate')
            nTrialsThisDate = [celldbThisDate.preToneNtrials.iloc[0],
                               celldbThisDate.salineToneNtrials.iloc[0],
                               celldbThisDate.doiToneNtrials.iloc[0]]
            nTrialsEachRunningCond[indc].append(nTrialsThisDate)
# -- Fraction of running trials under each reagent condition. --
fractionRunning = ( np.array(nTrialsEachRunningCond[1]) /
                    (np.array(nTrialsEachRunningCond[0])+np.array(nTrialsEachRunningCond[1])) )
pcRunning = 100*fractionRunning
modIndex1 = studyutils.modulation_index(pcRunning[:,1], pcRunning[:,0])
modIndex2 = studyutils.modulation_index(pcRunning[:,2], pcRunning[:,1])
#wstatMI, pValModInd = stats.wilcoxon(modIndex1, modIndex2)
wstatMI, pValModInd = (0,0)  # DEBUG
medianMI1 = np.nanmedian(modIndex1)
medianMI2 = np.nanmedian(modIndex2)
print(f'Running (%) MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')
print(f'Running (%) Median : pre={np.nanmedian(pcRunning[:,0]):0.1f},  ' +
      f'saline={np.nanmedian(pcRunning[:,1]):0.1f},  ' +
      f'DOI={np.nanmedian(pcRunning[:,2]):0.1f}')
print('')


# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
#gsMain = gridspec.GridSpec(1, len(studyparams.REAGENTS))
gsMain = gridspec.GridSpec(len(metrics), len(runningConditions)+1)
gsMain.update(left=0.05, right=0.99, top=0.92, bottom=0.08, wspace=0.3, hspace=0.4)
axs = []

# -- Show panel labels --
plt.figtext(labelPosX[0], labelPosY[0], 'A', fontsize=fontSizePanel, fontweight='bold')
plt.figtext(labelPosX[0], labelPosY[1], 'B', fontsize=fontSizePanel, fontweight='bold')
for indp, plabel in enumerate(['C','D','E']):
    plt.figtext(labelPosX[indp+1], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')

# -- Plot number of running trials for each reagent --
nSessions = len(pcRunning)
axTrials = plt.subplot(gsMain[1, 0])
plt.plot([0, 1, 2], pcRunning.T, '-', color='0.75')
#plt.plot([0, 1, 2], pcRunning.T, 'o', mec=pColor, color=pColor)
for indr, reagent in enumerate(studyparams.REAGENTS):
    thisColorLight = figparams.colorsLight[reagent]
    thisColorDark = figparams.colors[reagent]
    plt.plot(np.tile(indr, nSessions), pcRunning[:,indr], 'o', mec=thisColorLight, mfc=thisColorLight)
    medianFR = np.nanmedian(pcRunning[:,indr])
    plt.plot(indr, medianFR, '_', ms=40, mew=3, color=thisColorDark)
plt.xlim([-0.5, 2.5])
plt.ylim([0, 100])
plt.ylabel(f'Percent trials running (%)', fontsize=fontSizeLabels)
plt.xticks([0, 1, 2], studyparams.REAGENTS, fontsize=fontSizeTicks)
#wstat1, pVal1 = stats.wilcoxon(pcRunning[:,0], pcRunning[:,1])
#wstat2, pVal2 = stats.wilcoxon(pcRunning[:,1], pcRunning[:,2])
wstat1, pVal1 = (0,0)  # DEBUG
wstat2, pVal2 = (0,0)  # DEBUG
plt.text(0.5, 0.95, f'N = {nSessions} sessions\n'+f'({len(studyparams.SUBJECTS)} mice)',
         transform=axTrials.transAxes, ha='center', fontsize=fontSizeLabels)
extraplots.boxoff(axTrials)


# -- Plot baseline and evoked firing rates for each running condition --
for indrun, runningCond in enumerate(runningConditions):
    celldb = celldbs[indrun]
    for indm, metric in enumerate(metrics):

        thisMetricFiring = np.empty((len(celldb), len(studyparams.REAGENTS)))
        for indr, reagent in enumerate(studyparams.REAGENTS):
            thisMetricFiring[:, indr] = celldb[reagent + metric]

        #thisMetricFiring = thisMetricFiring[responsive & steady & posResponse, :] 
        thisMetricFiring = thisMetricFiring[responsive & steady, :]   # ORIGINAL
        #thisMetricFiring = thisMetricFiring[responsive, :]  # Control to compare pre-sal to sal-doi


        axs.append(plt.subplot(gsMain[indm, indrun+1]))
        nCells = thisMetricFiring.shape[0]
        thisMetricSaline = thisMetricFiring[:,1]
        thisMetricDOI = thisMetricFiring[:,2]
        #thisMetricDOI = thisMetricFiring[:,0]  # Pre
        notNan = ~(np.isnan(thisMetricSaline)|np.isnan(thisMetricDOI))
        wstat, pVal = stats.wilcoxon(thisMetricSaline[notNan],
                                     thisMetricDOI[notNan])
        medianSaline = np.nanmedian(thisMetricSaline)
        medianDOI = np.nanmedian(thisMetricDOI)
        print(f'[N={nCells}] {metric} Median : saline={medianSaline:0.4f}, ' +
              f'DOI={medianDOI:0.4f}   p={pVal:0.6f}')
        print(f'\t DOI>saline: {np.mean(thisMetricDOI>thisMetricSaline):0.1%}'+
              f'\t DOI<saline: {np.mean(thisMetricDOI<thisMetricSaline):0.1%}')
    
        plt.plot(thisMetricSaline, thisMetricDOI, 'o', mec='0.75', mfc='none')
        #maxVal = max(thisMetricSaline.max(), thisMetricDOI.max())
        maxVal = metricLims[indm][-1]
        plt.plot([0, maxVal], [0, maxVal], 'k-', lw=0.5)
        plt.plot(medianSaline, medianDOI, '+', color='m', ms=10, mew=2)
        plt.gca().set_aspect('equal', 'box')
        plt.xlim(metricLims[indm])
        plt.ylim(metricLims[indm])
        #plt.xlim([0, maxVal])
        #plt.ylim([0, maxVal])
        #plt.xticks(metricTicks[indm], fontsize=fontSizeTicks)
        #plt.yticks(metricTicks[indm], fontsize=fontSizeTicks)
        plt.xlabel(f'{metricLabels[indm]} Saline (spk/s)', fontsize=fontSizeLabels)
        plt.ylabel(f'{metricLabels[indm]} DOI (spk/s)', fontsize=fontSizeLabels)
        if indm == 0:
            plt.title(runningCondLabels[indrun]+f' ({nCells} cells)', fontsize=fontSizeLabels, fontweight='bold')

        #plt.text(metricLims[indm][-1]/2, 0.9*metricLims[indm][-1], f'p = {pVal:0.4f}',
        #         transform=axTrials.transAxes, ha='center', fontsize=fontSizeTicks)
        pValStr = f'p = {pVal:0.4f}' if pVal>=0.001 else 'p < 0.001'
        plt.text(0.5, 0.9, pValStr,
                 transform=axs[-1].transAxes, ha='center', fontsize=fontSizeTicks)
        #f'p = {pVal:0.4f}'
        modIndex1 = studyutils.modulation_index(thisMetricFiring[:,1], thisMetricFiring[:,0])
        modIndex2 = studyutils.modulation_index(thisMetricFiring[:,2], thisMetricFiring[:,1])
        hasModIndex = (~np.isnan(modIndex1) & ~np.isnan(modIndex2))
        wstatMI, pValModInd = stats.wilcoxon(modIndex1[hasModIndex], modIndex2[hasModIndex])
        medianMI1 = np.nanmedian(modIndex1)
        medianMI2 = np.nanmedian(modIndex2)
        print(f'[N={nCells}] {metric} MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')
        
        plt.show()
    #sys.exit()
    print('')
        

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
