"""
Figure about photoidentification.
"""

import sys
import studyparams
import studyutils
import figparams
import os
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from jaratoolbox import spikesorting
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from importlib import reload
reload(figparams)
reload(spikesanalysis)


SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'fig4_response_dynamics' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7, 2.8] # In inches

PANELS = [1, 1,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 1  #figparams.rasterMarkerSize

#labelPosX = [0.01, 0.22, 0.41, 0.59, 0.8]   # Horiz position for panel labels
labelPosX = [0.01, 0.27, 0.515, 0.78]   # Horiz position for panel labels
labelPosY = [0.94, 0.48]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
laserColor = figparams.colors['blueLaser']
soundColor = figparams.colors['soundStim']
stimLineWidth = 6

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning_latency.h5')
celldb = celldatabase.load_hdf(dbPath)

selcells = []

selcells.append(['d1pi041', '2019-08-27', 3400, 1, 3])  # cell 1547
selcells.append(['d1pi042', '2019-09-11', 3000, 4, 3])  # cell 1902
#selcells.append(['d1pi046', '2020-02-18', 3500, 7, 5])  # cell 3068
selcells.append(['d1pi041', '2019-08-25', 2700, 8, 2])  # cell 1188

selcells.append(['d1pi041', '2019-08-25', 3000, 7, 3])  # cell 1246
selcells.append(['d1pi042', '2019-09-11', 3400, 4, 5])  # cell 1999
selcells.append(['d1pi041', '2019-08-25', 2900, 5, 3])  # cell 1210

cellTypes = ['D1', 'D1', 'D1',  'ND1', 'ND1', 'ND1']
panelLabels = ['A', 'B', 'C', 'D', 'E', 'F']
axPos = [0,1,2, 3,4,5]  #5,6,7]
cellTypeLabel = ['D1', 'non-D1']
#selectedCells = [selectedD1, selectedND1]

indsToPlot = slice(0,12) # onset
indsToPlot = slice(9,21) # sustained
indsToPlot = slice(21,29) # delayed
indsToPlot = slice(14,26) # delayed

indsToPlot = slice(0,8) # delayed

#indsToPlot = slice(10,11)
#indsToPlot = slice(4,5)
selectedCells = selcells[indsToPlot]
cellTypes = cellTypes[indsToPlot]



# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(1, 2, width_ratios=[0.85, 0.15])
#gsMain.update(left=0.10, right=0.99, top=0.98, bottom=0.11, wspace=0.52, hspace=0.5)
gsMain.update(left=0.07, right=0.99, top=0.95, bottom=0.15, wspace=0.4, hspace=0.5)
gsRasters = gsMain[0,0].subgridspec(2, 3, wspace=0.4, hspace=0.25)
#gsBoxes = gsMain[0,1].subgridspec(1, 2, wspace=1)

def plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, color='k'):
    binEdges = np.arange(timeRange[0], timeRange[-1], 0.01)
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                             indexLimitsEachTrial, binEdges)
    binWidth = binEdges[1]-binEdges[0]
    firingRateMat = spikeCountMat/binWidth
    smoothPSTH = True
    smoothWinSize = 1
    extraplots.plot_psth(firingRateMat, smoothWinSize, binEdges, trialsEachCond=[],
                         colorEachCond=[color], linewidth=2, downsamplefactor=1)
    extraplots.boxoff(axPSTH)
    ylims = plt.ylim()
    yticks = [0, np.ceil(ylims[1])]
    plt.ylim(yticks)
    plt.yticks(yticks)
    #plt.ylabel('spk/s', fontsize=fontSizeLabels)

def plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, stimColor='g'):
    plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=markerSize)
    #plt.xlim(timeRange)
    plt.axis(False)
    yLims = plt.ylim()
    # -- Plot the stimulus --
    yPos = 1.1*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    plt.plot([0, 0.1], 2*[yPos], lw=stimLineWidth, color=stimColor,
             clip_on=False, solid_capstyle='butt')

# -- Load laser session --
if PANELS[0]:
    for indType, cellType in enumerate(cellTypes):
        colorThisType = figparams.colors[cellType]
        indRow, dbRow = celldatabase.find_cell(celldb, *selectedCells[indType])
        oneCell = ephyscore.Cell(dbRow)

        #gs00 = gsMain[indType].subgridspec(2, 1, hspace=0)

        # ---------------------- SOUND -------------------------
        freqIndex = dbRow.toneIndexBest
        #freqIndex = 8
        if 0:#indRow==3068:
            freqIndex = 8

        ephysData, bdata = oneCell.load('tuningCurve')
        spikeTimes = ephysData['spikeTimes']
        #eventOnsetTimes = ephysData['events']['stimOn']
        detectorOnsetTimes = ephysData['events']['soundDetectorOn']
        eventOnsetTimes = spikesanalysis.minimum_event_onset_diff(detectorOnsetTimes, 0.5)

        soundEventID = bdata.stateMatrix['statesNames']['output1On']
        soundEvents = bdata.events['nextState'] == soundEventID
        eventOnsetTimesBehav = bdata.events['eventTime'][soundEvents]
        missingTrials = behavioranalysis.find_missing_trials(eventOnsetTimes, eventOnsetTimesBehav)

        if len(bdata['currentFreq']) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(bdata['currentFreq'])]
        if len(bdata['currentFreq']) == len(eventOnsetTimes)-2:
            print(f'{indRow} was two trials off')
            ###bdata['currentFreq'] = bdata['currentFreq'][1:]
            offset = 1
            eventOnsetTimes = eventOnsetTimes[offset:len(bdata['currentFreq'])+offset]
            freqIndex = 8
        #timeRange = [-0.3, 0.6]  # In seconds
        timeRange = [-0.2, 0.4]  # In seconds

        thisFreq = np.unique(bdata['currentFreq'])[freqIndex]
        selectedTrials = (bdata['currentFreq']==thisFreq) & (bdata['currentIntensity']>40)

        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials],
                                                  timeRange)

        # -- Calculate latency --
        respLatency = dbRow.toneLatencySoundDetector

        # -- Plot raster/PSTH --
        gs01 = gsRasters[axPos[indType]].subgridspec(2, 1, hspace=0)
        axRaster = plt.subplot(gs01[0, 0])
        axPSTH = plt.subplot(gs01[1, 0], sharex=axRaster)
        labelRow = 0 if indType<3 else 1
        axRaster.annotate(panelLabels[indType], xy=(labelPosX[indType%3],labelPosY[labelRow]),
                          xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

        # -- Plot raster --
        plt.sca(axRaster)
        plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, soundColor)

        # -- Plot PSTH --
        plt.sca(axPSTH)
        plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, colorThisType)
        plt.xlim([-0.075, 0.25])
        plt.xticks(np.arange(-0.0, 0.3, 0.1))
        #plt.xlabel(f'{indRow} ({freqIndex}) Time from sound onset (s)', fontsize=fontSizeLabels)
        if indType in [3,4,5]:
            #plt.xlabel(f'Time from sound onset (s)', fontsize=fontSizeLabels)
            plt.xlabel(f'Time (s)', fontsize=fontSizeLabels)
        else:
            axPSTH.set_xticklabels([])
        if indType in [0, 3]:
            labelpad = 4 if plt.ylim()[1]<100 else -1
            plt.ylabel(f'spk/s', fontsize=fontSizeLabels, labelpad=labelpad)
        extraplots.set_ticks_fontsize(axPSTH, fontSizeTicks)
        #plt.axvline(0, color='k')
        #plt.axvline(respLatency, color='0.5')

# -- Select cells to compare --
restrictND1 = 0  # True: use only ND1 from tetrodes with D1
toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1 = \
    studyutils.select_cells(celldb, restrictND1=restrictND1)
cellsWithTone = ~np.isnan(celldb['toneMinPval'])
toneRespIndex = ( (celldb.toneFiringRateBest - celldb.toneFiringRateBaseline) /
                  (celldb.toneFiringRateBest + celldb.toneFiringRateBaseline) )
goodFit = celldb.toneGaussianRsquare > 0.01 #0.01

# -- Latency of selected cells --
#haveBW40 = ~np.isnan(celldb.toneBW40)
lowThresh = celldb.toneIntensityThreshold < 60  # (<55, l>0) (<60, l>0.005)
validLatency = celldb.toneLatencySoundDetector > 0.000
cellsToCompare = toneResponsive & validLatency & lowThresh

latencyD1 = celldb.toneLatencySoundDetector[indD1 & cellsToCompare]
latencyND1 = celldb.toneLatencySoundDetector[indND1 & cellsToCompare]
#latencyD1 = celldb.toneLatencySoundDetector[indD1 & toneResponsive & validLatency & haveBW40]
#latencyND1 = celldb.toneLatencySoundDetector[indND1 & toneResponsive & validLatency & haveBW40]

uval, pValLatency = stats.mannwhitneyu(latencyD1, latencyND1, alternative='two-sided')
print(f'Median tuning latency: D1 ({len(latencyD1)}): {1e3*np.median(latencyD1):0.2f} ms vs ' +
      f'ND1 ({len(latencyND1)}): {1e3*np.median(latencyND1):0.2f} ms   ' +
      f'p = {pValLatency:0.6f} (U={uval})')
print(f'Mean tuning latency: D1 ({len(latencyD1)}): {1e3*np.mean(latencyD1):0.2f} ms vs ' +
      f'ND1 ({len(latencyND1)}): {1e3*np.mean(latencyND1):0.2f} ms   ')

# -- Plot latency --
'''
axLatency = plt.subplot(gsBoxes[0])
axLatency.annotate('G', xy=(labelPosX[3],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
boxWidth = 0.5
lineWidth = 1.2
pboxD1 = plt.boxplot(1e3*latencyD1, positions=[1], widths=[boxWidth])
pboxND1 = plt.boxplot(1e3*latencyND1, positions=[2], widths=[boxWidth])
plt.setp(pboxD1['boxes'][0], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['medians'][0], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['whiskers'], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['caps'], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['fliers'], mec=figparams.colors['D1'], ms=3)
plt.setp(pboxND1['boxes'][0], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['medians'][0], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['whiskers'], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['caps'], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['fliers'], mec=figparams.colors['ND1'], ms=3)
plt.ylim([-0.2, 60])
plt.xticks([1,2])
plt.yticks(np.arange(0, 70, 10))
axLatency.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axLatency)
plt.ylabel('Response latency (ms)', fontsize=fontSizeLabels)
'''

# -------------------- Onset vs Sustain ------------------------
baselineFiring = celldb['toneFiringRateBaseline']
onsetResp = celldb['toneFiringRateBestOnset']
sustainResp = celldb['toneFiringRateBestSustain']
onsetToSustainIndex = (onsetResp - sustainResp) / (onsetResp + sustainResp)

#toneResponsive = toneResponsive & goodFit; print('WARNING! only good fit neurons used here!')
#toneResponsive = toneResponsive & ((celldb['toneFiringRateBest']-celldb['toneFiringRateBaseline'])>0)
#cellsToCompare = toneResponsive

# NOTE: I need to remove nans because some neurons may be have FR=0 for tones but not AM
osiD1 = onsetToSustainIndex[indD1 & cellsToCompare & ~np.isnan(onsetToSustainIndex)]
osiND1 = onsetToSustainIndex[indND1 & cellsToCompare & ~np.isnan(onsetToSustainIndex)]
nD1osi = len(osiD1)
nND1osi = len(osiND1)
medianD1osi = np.median(osiD1)
medianND1osi = np.median(osiND1)

uval, pValOSI = stats.mannwhitneyu(osiD1, osiND1, alternative='two-sided')
print(f'Onset-to-sustain index: D1 ({nD1osi}): {medianD1osi:0.4f} vs ' +
      f'ND1 ({nND1osi}): {medianND1osi:0.4f}   ' +
      f'p = {pValOSI:0.6f} (U={uval})\n')

# -- Plot Onset-to-Sustained ratio --
axOSI = plt.subplot(gsMain[1])
axOSI.annotate('G', xy=(labelPosX[3],labelPosY[0]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
boxWidth = 0.5
lineWidth = 1.2
#plt.axhline(0, color='0.5', ls='--')
pboxD1 = plt.boxplot(osiD1, positions=[1], widths=[boxWidth])
pboxND1 = plt.boxplot(osiND1, positions=[2], widths=[boxWidth])
plt.setp(pboxD1['boxes'][0], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['medians'][0], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['whiskers'], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['caps'], color=figparams.colors['D1'], lw=lineWidth)
plt.setp(pboxD1['fliers'], mec=figparams.colors['D1'], ms=3)
plt.setp(pboxND1['boxes'][0], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['medians'][0], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['whiskers'], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['caps'], color=figparams.colors['ND1'], lw=lineWidth)
plt.setp(pboxND1['fliers'], mec=figparams.colors['ND1'], ms=3)
plt.ylim([-1, 1])
plt.xticks([1,2])
plt.yticks(np.arange(-1, 1.5, 0.5))
axOSI.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axOSI)
plt.ylabel('Onset-to-sustained index', fontsize=fontSizeLabels, labelpad=3)


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)


sys.exit()


# -- Find low latency cells --
selD1 = (indD1 & toneResponsive & validLatency & lowThresh &
         (celldb.toneLatencySoundDetector>0) & (celldb.toneLatencySoundDetector<0.01)) 
selD1x = (indD1 & toneResponsive & validLatency & lowThresh &
         (celldb.toneLatencySoundDetector>0.04)) 
celldb[selD1].toneLatencySoundDetector*1e3
selND1 = (indND1 & toneResponsive & validLatency & lowThresh &
         (celldb.toneLatencySoundDetector>0) & (celldb.toneLatencySoundDetector<0.01)) 
selND1x = (indND1 & toneResponsive & validLatency & lowThresh &
         (celldb.toneLatencySoundDetector>0.04)) 
celldb[selND1].toneLatencySoundDetector*1e3
