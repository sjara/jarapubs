"""
This script plots changes in freq tuning from saline to DOI.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
from jaratoolbox import behavioranalysis
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_freqtuning' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [14, 4.8] # In inches

if len(sys.argv)==2:
    trialSubset = sys.argv[1]
    figFilename = figFilename + '_' + trialSubset
    figFormat = 'pdf'
else:
    trialSubset = ''
if trialSubset not in ['', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

# -- Example cells to plot --
#Cells with decrease: 149, 200, 430, 432, 669, 732, 806, 1145 
cellsToPlot = [149, 125, 83]   # Best: 149, 200, 1145
#cellsToPlot = [1145, 125, 83]   # Best: 149, 200, 1145
#cellsToPlot = [407, 1320, 1474]  # Give negative estimated max firing rate on DOI

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.002, 0.215, 0.43, 0.65, 0.825]   # Horiz position for panel labels
labelPosY = [0.96, 0.47]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorsRasterDark = figparams.colors
colorsRasterLight = figparams.colorsLight
pColor = '0.5'
#colorsRasterDark = {'saline': figparams.colors['saline'], 'doi': figparams.colors['doi']} 
#colorsRasterLight = {'saline': figparams.colorsLight['saline'], 'doi': figparams.colorsLight['doi']} 
#Oddball = figparams.colors['oddball']
#colorStandard = figparams.colors['standard']

rasterMarkerSize = 0.5

def plot_stim(yLims, stimDuration, stimLineWidth=6, stimColor=figparams.colorStim):
    yPos = 1.02*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    pstim = plt.plot([0, stimDuration], 2*[yPos], lw=stimLineWidth, color=stimColor,
                     clip_on=False, solid_capstyle='butt')
    return pstim[0]

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldbAll = celldatabase.load_hdf(dbFilename)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
    celldb = celldatabase.load_hdf(dbFilename)
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
    celldb = celldatabase.load_hdf(dbFilename)
else:
    celldb = celldbAll

# -- Process data --
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR

#responsive = studyutils.find_tone_responsive_cells(celldb, frThreshold=5)
#selective = studyutils.find_freq_selective(celldb, minR2=studyparams.MIN_R_SQUARED)
goodFit = studyutils.find_good_gaussian_fit(celldb, minR2=studyparams.MIN_R_SQUARED)
anyFitDOI = ~np.isnan(celldb['doiToneGaussianA'])
for indr, reagent in enumerate(studyparams.REAGENTS):
    #celldb[reagent+'ToneGaussianMax'] = ( celldb[reagent+'ToneGaussianA'] +
    #                                      celldb[reagent+'ToneGaussianY0'] )
    negResponseThisReagent = celldb[reagent+'ToneGaussianA']<0
    thisToneGaussianMax = ( celldb[reagent+'ToneGaussianA'] +
                            celldb[reagent+'ToneGaussianY0'] )
    thisToneGaussianMax[negResponseThisReagent] = celldb[reagent+'ToneGaussianY0']
    celldb[reagent+'ToneGaussianMax'] = thisToneGaussianMax
    celldb[reagent+'ToneGaussianBandwidth'] = \
        extraplots.gaussian_full_width_half_max(celldb[reagent+'ToneGaussianSigma'])
    baselineFiringRate = celldb[reagent+'ToneBaselineFiringRate']
    celldb[reagent+'ToneGaussianMaxChange'] = np.abs(thisToneGaussianMax-baselineFiringRate)
#posResponse = (celldbAll['preToneGaussianA']>0) #& (celldbAll['salineToneGaussianA']>0) #& (celldbAll['doiToneGaussianA']>0)
#posResponse = (celldbAll['salineToneGaussianA']>0) #& (celldbAll['doiToneGaussianA']>0)
#posResponse = (celldbAll['doiToneGaussianA']>0)
#posResponse = (celldb['salineToneGaussianMax']>0) & (celldb['doiToneGaussianMax']>0)
posResponse = (celldbAll['preToneGaussianA']>0)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldb, steadyParams, maxChangeFactor)

'''
metrics = ['ToneGaussianX0', 'ToneGaussianBandwidth', 'ToneGaussianMax']
metricsLabel = ['BF', 'Bandwidth', 'Max resp']
metricsUnits = ['kHz', 'oct', 'spk/s']
metricLims = [[np.log2(2000), np.log2(40000)], [0, 10], [0, 120]]
metricTicks = [[2000, 8000, 32000], [0, 5, 10], [0, 60, 120]]
'''

metrics = ['ToneGaussianX0', 'ToneGaussianBandwidth', 'ToneGaussianMaxChange']
metricsLabel = ['BF', 'Width', 'Max Δ']
#metricsLabel = ['BF', 'Bandwidth', 'Resp at BF']
metricsUnits = ['kHz', 'norm.', 'spk/s']
metricLims = [[np.log2(2000), np.log2(40000)], [0, 10], [0, 100]]
metricTicks = [[2000, 8000, 32000], [0, 5, 10], [0, 50, 100]]

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
gsMain = gridspec.GridSpec(1, 4, width_ratios=[0.2, 0.2, 0.2, 0.4])
gsMain.update(left=0.04, right=0.99, top=0.95, bottom=0.1, wspace=0.35, hspace=0.3)

# -- Show panel labels --
for indp, plabel in enumerate(['A','B','C','D','E']):
    plt.figtext(labelPosX[indp], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')
for indp, plabel in enumerate(['F','G']):
    plt.figtext(labelPosX[indp+3], labelPosY[1], plabel, fontsize=fontSizePanel, fontweight='bold')

# -- Raster and PSTH parameters --
timeRange = [-0.2, 0.4]
binWidth = 0.010
timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
smoothWinSizePsth = 2 
lwPsth = 2.5
downsampleFactorPsth = 1

# -- Plot examples --
if 1:
    sessionType = 'PureTones'
    for indcell, cellInd in enumerate(cellsToPlot):
        gsExample = gsMain[0, indcell].subgridspec(2, 1, hspace=0.35)
        gsRasters = gsExample[0].subgridspec(2, 1, hspace=0.1)
        axTuning = plt.subplot(gsExample[1])
        #cellInd, dbRow = celldatabase.find_cell(celldb, **cellsToPlot[indcell])
        dbRow = celldb.iloc[cellsToPlot[indcell]]
        oneCell = ephyscore.Cell(dbRow)
        reagentsToPlot = ['saline', 'doi']
        reagentLabels = ['Saline', 'DOI']
        allFits = []
        
        for indr, reagent in enumerate(reagentsToPlot):
            ephysData, bdata = oneCell.load(reagent+sessionType)  
            spikeTimes = ephysData['spikeTimes']
            eventOnsetTimes = ephysData['events']['stimOn']

            stimEachTrial = bdata['currentFreq']
            nTrials = len(stimEachTrial)
            
            # -- Doing this before selecting trials by running/notrunning --
            possibleStim = np.unique(stimEachTrial)
            nStim = len(possibleStim)
            
            stimDuration = bdata['stimDur'][-1]

            # If the ephys data is 1 more than the bdata, delete the last ephys trial.
            if len(stimEachTrial) == len(eventOnsetTimes)-1:
                eventOnsetTimes = eventOnsetTimes[:len(stimEachTrial)]
            assert len(stimEachTrial) == len(eventOnsetTimes), \
                "Number of trials in behavior and ephys do not match for {oneCell}"

            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

            trialsEachCond = behavioranalysis.find_trials_each_type(stimEachTrial, possibleStim)

            # -- Estimate evoked firing rate for each stim --
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial,
                                                                     timeRange)

            # -- Plot Raster --
            axRaster = plt.subplot(gsRasters[indr])
            possibleStimInKHz = possibleStim/1000
            rasterLabels = ['']*nStim;
            rasterLabels[0] = int(possibleStimInKHz[0])
            rasterLabels[-1] = int(possibleStimInKHz[-1])
            colorEachCond = [colorsRasterDark[reagent], colorsRasterLight[reagent]]*(nStim//2+1)
            (pRasterS,hcond,zline) = extraplots.raster_plot(spikeTimesFromEventOnset,
                                                            indexLimitsEachTrial, timeRange,
                                                            trialsEachCond, labels=rasterLabels,
                                                            colorEachCond=colorEachCond,
                                                            rasterized=True)
            plt.setp(pRasterS, ms=rasterMarkerSize)
            plt.xlabel('Time (s)', fontsize=fontSizeLabels)
            plt.ylabel('Freq (kHz)', fontsize=fontSizeLabels)
            #axRaster.set_yticklabels(['2']+['']*(nFreq-2)+['40'])
            if indr==0:
                plot_stim(plt.ylim(), stimDuration)
                axRaster.set_xticklabels([])

            # -- Plot tuning curve --
            plt.sca(axTuning)
            firingRates = dbRow[reagent+'ToneFiringRateEachFreq']
            fitParams = [dbRow[reagent+'ToneGaussianA'], dbRow[reagent+'ToneGaussianX0'],
                         dbRow[reagent+'ToneGaussianSigma'], dbRow[reagent+'ToneGaussianY0']]
            pdots, pfit = extraplots.plot_tuning_curve(possibleStim, firingRates, fitParams)
            pfit[0].set_color(colorsRasterDark[reagent])
            allFits.append(pfit[0])
            pdots[0].set_color(colorsRasterDark[reagent])
            extraplots.boxoff(axTuning)
            xTicks = np.array([2000, 4000, 8000,16000, 32000])
            axTuning.set_xticks(np.log2(xTicks))
            axTuning.set_xticklabels((xTicks/1000).astype(int))
            axTuning.set_xlabel('Frequency (kHz)', fontsize=fontSizeLabels)
            axTuning.set_ylabel('Firing rate (spk/s)', fontsize=fontSizeLabels)
        axTuning.legend([allFits[0], allFits[1]], ['Saline', 'DOI'], loc='upper left', handlelength=1)


# -- Plot metrics summary --        
gsMetrics = gsMain[0, 3].subgridspec(2, 2, hspace=0.4, wspace=0.5)
#for indm, metric in enumerate(['ToneGaussianMax']): #enumerate(metrics):
for indm, metric in enumerate(metrics):
    selectedCells = steady & goodFit & posResponse & anyFitDOI
    nCells = np.sum(selectedCells)
    thisMetricSaline = celldb['saline'+metric][selectedCells]
    thisMetricDOI = celldb['doi'+metric][selectedCells]
    wstat, pVal = stats.wilcoxon(thisMetricSaline, thisMetricDOI)
    medianSaline = np.median(thisMetricSaline)
    medianDOI = np.median(thisMetricDOI)
    print(f'[N={nCells}]  Median {metric}: saline={medianSaline:0.4f}, ' +
          f'DOI={medianDOI:0.4f}   p={pVal:0.4f}')
    print(f'\t DOI>saline: {np.mean(thisMetricDOI>thisMetricSaline):0.1%}'+
          f'\t DOI<saline: {np.mean(thisMetricDOI<thisMetricSaline):0.1%}')
    
    thisAx = plt.subplot(gsMetrics[indm+1])
    maxVal = max(thisMetricSaline.max(), thisMetricDOI.max())
    if indm==99:
        plt.loglog(thisMetricSaline, thisMetricDOI, 'o', mec='0.75', mfc='none')
    else:
        plt.plot(thisMetricSaline, thisMetricDOI, 'o', mec='0.75', mfc='none')
    plt.plot([0, maxVal], [0, maxVal], 'k-', lw=0.5)
    plt.plot(medianSaline, medianDOI, '+', color='m', ms=10, mew=2)
    plt.gca().set_aspect('equal', 'box')
    plt.xlim(metricLims[indm])
    plt.ylim(metricLims[indm])
    if indm==0:
        axTicksLog = np.log2(metricTicks[indm])
        axTickLabels = [f'{int(x/1000)}' for x in metricTicks[indm]]
        plt.xticks(axTicksLog, axTickLabels, fontsize=fontSizeTicks)
        plt.yticks(axTicksLog, axTickLabels, fontsize=fontSizeTicks)
    else:
        plt.xticks(metricTicks[indm], fontsize=fontSizeTicks)
        plt.yticks(metricTicks[indm], fontsize=fontSizeTicks)
    plt.xlabel(f'{metricsLabel[indm]} Saline ({metricsUnits[indm]})', fontsize=fontSizeLabels)
    plt.ylabel(f'{metricsLabel[indm]} DOI ({metricsUnits[indm]})', fontsize=fontSizeLabels)
    #plt.gca().tick_params(labelleft=False)
    if indm==99:
        thisAx.set_xscale('log')
        thisAx.set_yscale('log')
    plt.text(0.5, 0.9, f'p = {pVal:0.3f}',
             transform=thisAx.transAxes, ha='center', fontsize=fontSizeTicks)
    if indm==0:
        plt.title(f'N = {nCells} cells', fontsize=fontSizeLabels, fontweight='normal')
    #extraplots.boxoff(thisAx)
    plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
