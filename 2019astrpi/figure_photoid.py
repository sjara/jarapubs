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


SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_photoid' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
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
dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220129.h5')
celldb = celldatabase.load_hdf(dbPath)


selectedD1 = ['d1pi041', '2019-08-27',3400, 4, 6]  # cell 1557
selectedND1 = ['d1pi041', '2019-08-25', 3500, 7, 3] # cell 1421
cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
selectedCells = [selectedD1, selectedND1]

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 3)
gsMain.update(left=0.10, right=0.99, top=0.98, bottom=0.11, wspace=0.52, hspace=0.5)

axCartoon = plt.subplot(gsMain[0, 0])
axCartoon.set_axis_off()
axCartoon.annotate(panelLabels[0], xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')

def plot_spikeshapes(waveforms):
    (nSpikes, nChannels, nSamplesPerSpike) = waveforms.shape
    meanWaveforms = waveforms.mean(axis=0)
    stdWaveforms = waveforms.std(axis=0)
    xRange = np.arange(nSamplesPerSpike)
    for indc in range(nChannels):
        newXrange = xRange+indc*(nSamplesPerSpike+12)
        pMean = plt.plot(newXrange, meanWaveforms[indc,:], color='k', lw=2, clip_on=False)
        pStd = plt.fill_between(newXrange, meanWaveforms[indc,:]+stdWaveforms[indc,:],
                                     meanWaveforms[indc,:]-stdWaveforms[indc,:], color='0.8',
                                     clip_on=False)
    plt.axis(False)
    
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
    yPos = 1.1*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    plt.plot([0, 0.1], 2*[yPos], lw=stimLineWidth, color=stimColor,
             clip_on=False, solid_capstyle='butt')

# -- Load laser session --
for indType, cellType in enumerate(cellTypes):
    colorThisType = figparams.colors[cellType]
    indRow, dbRow = celldatabase.find_cell(celldb, *selectedCells[indType])
    oneCell = ephyscore.Cell(dbRow)
    
    # ---------------------- LASER -------------------------
    ephysData, bdata = oneCell.load('laserpulse')
    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    timeRange = [-0.3, 0.6]  # In seconds
    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    # -- Plot laser spikeShapes/raster/PSTH --
    gs00 = gsMain[0, indType+1].subgridspec(3, 1, hspace=0)
    axShapes = plt.subplot(gs00[0, 0])
    axRaster = plt.subplot(gs00[1, 0])
    axPSTH = plt.subplot(gs00[2, 0], sharex=axRaster)
    axShapes.annotate(panelLabels[indType+1][0], xy=(labelPosX[indType+1],labelPosY[0]),
                      xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

    # -- Plot waveforms --
    plt.sca(axShapes)
    plot_spikeshapes(ephysData['samples'])
    
    # -- Plot raster --
    plt.sca(axRaster)
    plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, laserColor)
    
    # -- Plot PSTH --
    plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, colorThisType)
    plt.xlim([-0.25, 0.45])
    plt.xticks([-0.2, 0, 0.2, 0.4])
    plt.ylim([0, 65])
    plt.text(0.3, 25, cellTypeLabel[indType], ha='center', fontweight='bold', color=colorThisType)
    plt.xlabel('Time from laser onset (s)', fontsize=fontSizeLabels)

    # ---------------------- SOUND -------------------------
    freqIndex = 8 # 9883.54 khz
    ephysData, bdata = oneCell.load('tuningCurve')
    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    if len(bdata['currentFreq']) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(bdata['currentFreq'])]
    timeRange = [-0.3, 0.6]  # In seconds

    thisFreq = np.unique(bdata['currentFreq'])[freqIndex]
    selectedTrials = (bdata['currentFreq']==thisFreq) & (bdata['currentIntensity']>50)

    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials],
                                              timeRange)
    
    # -- Plot laser spikeShapes/raster/PSTH --
    gs01 = gsMain[1, indType+1].subgridspec(3, 1, hspace=0)
    axShapes = plt.subplot(gs01[0, 0])
    axRaster = plt.subplot(gs01[1, 0])
    axPSTH = plt.subplot(gs01[2, 0], sharex=axRaster)
    axShapes.annotate(panelLabels[indType+1][1], xy=(labelPosX[indType+1],labelPosY[1]),
                      xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
    
    # -- Plot waveforms --
    plt.sca(axShapes)
    plot_spikeshapes(ephysData['samples']) #[trialsThisFreq, :, :]
    
    # -- Plot raster --
    plt.sca(axRaster)
    plot_raster(spikeTimesFromEventOnset, trialIndexForEachSpike, soundColor)
    
    # -- Plot PSTH --
    plot_psth(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, colorThisType)
    plt.xlim([-0.25, 0.45])
    plt.ylim([0, 36])
    plt.text(0.3, 15, cellTypeLabel[indType], ha='center', fontweight='bold', color=colorThisType)
    plt.xticks([-0.2, 0, 0.2, 0.4])
    plt.xlabel('Time from sound onset (s)', fontsize=fontSizeLabels)
    
plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
