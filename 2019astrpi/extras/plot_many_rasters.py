"""
Figure about photoidentification.
"""

import sys
import studyparams
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

'''
FIGNAME = 'figure_name'
figDataFile = 'file_containing_data_for_this_fig.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
'''

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_response_dynamics' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7,4] # In inches


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

panelLabels = ['A', ['B','D'], ['C','E']]
#labelPosX = [0.04, 0.3, 0.68]   # Horiz position for panel labels
labelPosX = [0.01, 0.37, 0.705]   # Horiz position for panel labels
labelPosY = [0.93, 0.41]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
laserColor = figparams.colors['blueLaser']
soundColor = figparams.colors['soundStim']
stimLineWidth = 6

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)

selcells = []

selcells.append(['d1pi032', '2019-02-19', 3200, 8, 4])  # cell 516
selcells.append(['d1pi041', '2019-08-27', 3400, 1, 3])  # cell 1547 * 
selcells.append(['d1pi041', '2019-08-27', 3400, 4, 6])  # cell 1557
selcells.append(['d1pi042', '2019-09-11', 3100, 4, 3])  # cell 1927 **
selcells.append(['d1pi042', '2019-09-11', 3200, 2, 2])  # cell 1944

selcells.append(['d1pi041', '2019-08-25', 3000, 7, 3])  # cell 1246 *
selcells.append(['d1pi041', '2019-08-25', 3100, 7, 3])  # cell 1278
selcells.append(['d1pi041', '2019-08-27', 3400, 4, 3])  # cell 1554
selcells.append(['d1pi049', '2020-03-13', 3200, 7, 3])  # cell 4416 **

selcells.append(['d1pi042', '2019-09-11', 3000, 4, 3])  # cell 1902 *
selcells.append(['d1pi042', '2019-09-11', 3200, 3, 4])  # cell 1949
selcells.append(['d1pi047', '2020-02-26', 3700, 4, 3])  # cell 3436
selcells.append(['d1pi036', '2019-05-29', 2800, 2, 4])  # cell 710
selcells.append(['d1pi036', '2019-05-29', 2900, 6, 4])  # cell 739

selcells.append(['d1pi042', '2019-09-11', 3400, 4, 5])  # cell 1999 *
selcells.append(['d1pi041', '2019-08-25', 3400, 7, 6])  # cell 1391
selcells.append(['d1pi041', '2019-08-25', 3600, 7, 2])  # cell 1454
selcells.append(['d1pi033', '2019-04-17', 2900, 8, 2])  # cell 650
selcells.append(['d1pi036', '2019-05-29', 2900, 5, 6])  # cell 736
selcells.append(['d1pi041', '2019-08-25', 3200, 8, 3])  # cell 1317
selcells.append(['d1pi041', '2019-08-25', 3300, 8, 6])  # cell 1359

selcells.append(['d1pi046', '2020-02-18', 3200, 6, 3])  # cell 2994
selcells.append(['d1pi046', '2020-02-18', 3300, 8, 5])  # cell 3021 *
selcells.append(['d1pi046', '2020-02-18', 3500, 7, 5])  # cell 3068 **
selcells.append(['d1pi046', '2020-02-18', 3700, 7, 2])  # cell 3126

selcells.append(['d1pi032', '2019-02-22', 2900, 8, 5])  # cell 537
selcells.append(['d1pi041', '2019-08-25', 2900, 5, 3])  # cell 1210 *
selcells.append(['d1pi043', '2020-01-15', 3200, 8, 6])  # cell 2101
selcells.append(['d1pi046', '2020-02-18', 3500, 8, 4])  # cell 3073

#selectedD1 = ['d1pi041', '2019-08-27',3400, 4, 6]  # cell 1557
#selectedND1 = ['d1pi041', '2019-08-25', 3500, 7, 3] # cell 1421
cellTypes = ['D1', 'D1', 'D1', 'D1', 'D1',  'ND1', 'ND1', 'ND1', 'ND1',
             'D1', 'D1', 'D1', 'D1', 'D1',  'ND1', 'ND1', 'ND1', 'ND1', 'ND1', 'ND1', 'ND1',
             'D1', 'D1', 'D1', 'D1',  'ND1', 'ND1', 'ND1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
#selectedCells = [selectedD1, selectedND1]

indsToPlot = slice(0,12) # onset
indsToPlot = slice(9,21) # sustained
indsToPlot = slice(21,29) # delayed
indsToPlot = slice(14,26) # delayed

#indsToPlot = slice(10,11)
#indsToPlot = slice(4,5)
selectedCells = selcells[indsToPlot]
cellTypes = cellTypes[indsToPlot]



# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(3, 4)
gsMain.update(left=0.10, right=0.99, top=0.98, bottom=0.11, wspace=0.52, hspace=0.5)

'''
axCartoon = plt.subplot(gsMain[0, 0])
axCartoon.set_axis_off()
axCartoon.annotate(panelLabels[0], xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
'''

def plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, color='k'):
    binEdges = np.arange(timeRange[0], timeRange[-1], 0.01)
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                             indexLimitsEachTrial, binEdges)
    binWidth = binEdges[1]-binEdges[0]
    firingRateMat = spikeCountMat/binWidth
    smoothPSTH = True
    smoothWinSize = 1
    plt.sca(axPSTH)
    extraplots.plot_psth(firingRateMat, smoothWinSize, binEdges, trialsEachCond=[],
                         colorEachCond=[color], linewidth=2.5, downsamplefactor=1)
    extraplots.boxoff(axPSTH)
    plt.ylabel('spk/s', fontsize=fontSizeLabels)

def plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, stimColor='g'):
    plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=markerSize)
    #plt.xlim(timeRange)
    plt.axis(False)
    yLims = plt.ylim()
    # -- Plot the stimulus --
    yPos = yLims[-1] + 0.1*(yLims[-1]-yLims[0])
    plt.plot([0, 0.1], 2*[yPos], lw=stimLineWidth, color=stimColor, clip_on=False)

# -- Load laser session --
for indType, cellType in enumerate(cellTypes):
    colorThisType = figparams.colors[cellType]
    indRow, dbRow = celldatabase.find_cell(celldb, *selectedCells[indType])
    oneCell = ephyscore.Cell(dbRow)
    
    #gs00 = gsMain[indType].subgridspec(2, 1, hspace=0)

    # ---------------------- SOUND -------------------------
    freqIndex = dbRow.toneIndexBest
    #freqIndex = 8
    
    ephysData, bdata = oneCell.load('tuningCurve')
    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']

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
    
    # -- Plot raster/PSTH --
    gs01 = gsMain[indType].subgridspec(2, 1, hspace=0)
    axRaster = plt.subplot(gs01[0, 0])
    axPSTH = plt.subplot(gs01[1, 0], sharex=axRaster)
    '''
    axRaster.annotate(panelLabels[indType+1][1], xy=(labelPosX[indType+1],labelPosY[1]),
                      xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
    '''
    
    # -- Plot raster --
    plt.sca(axRaster)
    plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, soundColor)
    
    # -- Plot PSTH --
    plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, colorThisType)
    #plt.xlim([-0.25, 0.45])
    plt.xlim([-0.2, 0.4])
    #plt.ylim([0, 36])
    #plt.text(0.3, 15, cellTypeLabel[indType], ha='center', fontweight='bold', color=colorThisType)
    plt.xticks([-0.2, 0, 0.2, 0.4])
    plt.xlabel(f'{indRow} ({freqIndex}) Time from sound onset (s)', fontsize=fontSizeLabels)
    
plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
