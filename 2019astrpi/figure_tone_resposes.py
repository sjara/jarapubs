"""
Figure about responses to pure tones.
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
import studyutils

from importlib import reload
reload(figparams)
reload(studyutils)

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_tone_responses' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [4.6, 3.8] # In inches
figSize = [7, 3.6] # In inches


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

#panelLabels = ['A', ['B','D'], ['C','E']]
#labelPosX = [0.04, 0.3, 0.68]   # Horiz position for panel labels
labelPosXtop = [0.01, 0.36, 0.68]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.55, 0.75]   # Horiz position for panel labels
labelPosY = [0.93, 0.5]    # Vert position for panel labels

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220212.h5')
dbPath = os.path.join(figuresDataDir, 'astrpi_am_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)

pp = lambda x: f"['{x.subject}', '{x.date}', {x.pdepth}, {x.egroup}, {x.cluster}]"
# pp(celldb.loc[])

# D1: 1967, 1971
selectedD1 = ['d1pi036', '2019-05-29', 2800, 2, 4] # cell 710
selectedND1 = ['d1pi041', '2019-08-25', 3600, 7, 2] # cell 1454
selectedND1sup = ['d1pi041', '2019-08-25', 3400, 6, 2] # cell 1381
#selectedND1sup = ['d1pi041', '2019-08-25', 3600, 5, 4] # cell 1446
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 2, 3]  # cell 1971
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 1, 2]  # cell 1967
#selectedND1 = ['d1pi041', '2019-08-25', 3500, 8, 2] # cell 1426
cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
selectedCells = [selectedD1, selectedND1, selectedND1sup]
examplesTuning = ['D1', 'ND1', 'ND1']
cMaps = ['Blues', 'Reds', 'Reds']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 1)
gsMain.update(left=0.1, right=0.95, top=1.0, bottom=0.12, hspace=0.3)

gsTC = gsMain[0,0].subgridspec(1, 3, wspace=0.6)
axTCD1 = plt.subplot(gsTC[0, 0])
axTCD1.annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCND1 = plt.subplot(gsTC[0, 1])
axTCND1.annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
axTCSup = plt.subplot(gsTC[0, 2])
axTCSup.annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
axTC = [axTCD1, axTCND1, axTCSup]

gsBottom = gsMain[1,0].subgridspec(1, 3, wspace=0.6, width_ratios=[0.6, 0.15, 0.15])
axHist = plt.subplot(gsBottom[0, 0])
axHist.annotate('D', xy=(labelPosXbot[0],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axBW = plt.subplot(gsBottom[0, 1])
axBW.annotate('E', xy=(labelPosXbot[1],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axOSI = plt.subplot(gsBottom[0, 2])
axOSI.annotate('F', xy=(labelPosXbot[2],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')


if 1:
    for indType, cellType in enumerate(examplesTuning):
        plt.sca(axTC[indType])
        #colorThisType = figparams.colors[cellType]
        indRow, dbRow = celldatabase.find_cell(celldb, *selectedCells[indType])
        oneCell = ephyscore.Cell(dbRow)

        # ---------------------- SOUND -------------------------
        ephysData, bdata = oneCell.load('tuningCurve')
        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        if len(bdata['currentFreq']) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(bdata['currentFreq'])]
        timeRange = [-0.3, 0.6]  # In seconds
        responsePeriod = [0, 0.1]
        respPeriodDuration = responsePeriod[1]-responsePeriod[0]

        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 responsePeriod)
        nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

        freqEachTrial = bdata['currentFreq']
        intensityEachTrial = bdata['currentIntensity']
        possibleFreq = np.unique(freqEachTrial)
        nFreq = len(possibleFreq)
        possibleIntensity = np.unique(intensityEachTrial)
        nIntensity = len(possibleIntensity)
        trialsEachComb = behavioranalysis.find_trials_each_combination(freqEachTrial,
                                                                       possibleFreq,
                                                                       intensityEachTrial,
                                                                       possibleIntensity)
        tuningMatSpikeCount = np.empty((nIntensity, nFreq))
        for indf, thisFreq in enumerate(possibleFreq):
            for indi, thisIntensity in enumerate(possibleIntensity):
                theseTrials = trialsEachComb[:,indf, indi]
                tuningMatSpikeCount[indi, indf] = np.mean(nSpikesEachTrial[theseTrials])
        tuningMatFiringRate = tuningMatSpikeCount/respPeriodDuration


        maxFR = np.max(tuningMatFiringRate)
        if indType==2:
            maxFR = 24; print('WARNING! fixing max color for third example.\n')
        tuningPlot = plt.imshow(np.flipud(tuningMatFiringRate), interpolation='nearest',
                                cmap=cMaps[indType], vmax=maxFR) #, extent=[0, 15, 11, 0]
        tuningPlot.set_clim([0, maxFR])

        #nIntenLabels = 3
        #intensities = np.linspace(possibleIntensity[0], possibleIntensity[-1], nIntenLabels)
        intensities = np.array([15, 42, 70])
        nIntenLabels = len(intensities)
        intenTickFactor = (nIntensity-1)/(intensities[-1]-intensities[0])
        intenTickLocations = intenTickFactor * (intensities[-1] - intensities)
        plt.yticks(intenTickLocations)
        axTC[indType].set_yticklabels(intensities)
        
        nFreqLabels = 3
        freqTickLocations = np.linspace(0, nFreq-1, nFreqLabels)
        #freqTickLocations = np.array([2, 20, 40])
        #nFreqLabels = len(freqTickLocations)
        #freqTickFactor = (nFreq-1)/(intensities[-1]-intensities[0])
        freqs = np.logspace(np.log10(possibleFreq[0]),np.log10(possibleFreq[-1]),nFreqLabels)
        #freqs = np.round(freqs, decimals=1)
        #freqLabels = ['{0:.1f}'.format(freq/1000) for freq in freqs]
        freqLabels = ['{0:0.0f}'.format(freq/1000) for freq in freqs]

        plt.xticks(freqTickLocations)
        if 1: #indType==1:
            axTC[indType].set_xticklabels(freqLabels)
            plt.xlabel('Frequency (kHz)', fontsize=fontSizeLabels)
        else:
            axTC[indType].set_xticklabels([])
        if indType==0:
            plt.ylabel('Intensity\n(dB SPL)', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(axTC[indType], fontSizeTicks)

        cbar = plt.colorbar(tuningPlot, ax=axTC[indType], format='%d', shrink=0.62,
                            pad=0.06, aspect=15)
        cbar.ax.set_ylabel('spk/s', fontsize=fontSizeLabels, labelpad=-2, rotation=-90)
        extraplots.set_ticks_fontsize(cbar.ax, fontSizeTicks)
        cbar.set_ticks([0, maxFR])

        plt.show()


# -- Tone response index --
restrictND1 = 0  # True: use only ND1 from tetrodes with D1
toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1 = \
    studyutils.select_cells(celldb, restrictND1=restrictND1)

cellsWithTone = ~np.isnan(celldb['toneMinPval'])
toneRespIndex = ( (celldb.toneFiringRateBest - celldb.toneFiringRateBaseline) /
                  (celldb.toneFiringRateBest + celldb.toneFiringRateBaseline) )
if 0:
    print('WARNING! only tone-responsive neurons used here!')
    goodFit = celldb.toneGaussianRsquare > 0.01
    #toneResponsive = toneResponsive & goodFit; print('WARNING! only good fit neurons used here!')
    toneRespIndexD1 = toneRespIndex[indD1 & toneResponsive]
    toneRespIndexND1 = toneRespIndex[indND1 & toneResponsive]
    toneRespIndexEachType = [toneRespIndex[indD1 & toneResponsive],
                             toneRespIndex[indND1 & toneResponsive]]
else:   
    toneRespIndexD1 = toneRespIndex[indD1 & cellsWithTone]
    toneRespIndexND1 = toneRespIndex[indND1 & cellsWithTone]
    toneRespIndexEachType = [toneRespIndexD1, toneRespIndexND1]

toneRespIndexD1pos = toneRespIndexD1[toneRespIndexD1>0]
toneRespIndexND1pos = toneRespIndexND1[toneRespIndexND1>0]
toneRespIndexD1neg = toneRespIndexD1[toneRespIndexD1<0]
toneRespIndexND1neg = toneRespIndexND1[toneRespIndexND1<0]
toneRespIndexEachTypePos = [toneRespIndexD1pos, toneRespIndexND1pos]
toneRespIndexEachTypeNeg = [toneRespIndexD1neg, toneRespIndexND1neg]


def line_histogram(data, edges, percent=False, **kwargs):
    hh, ee = np.histogram(data, edges)
    if percent:
        hh = 100 * hh/len(data[~np.isnan(data)])
    #print(np.sum(hh))
    hline = plt.step(np.r_[ee, ee[-1]], np.r_[0,hh,0],'-', where='pre', **kwargs)
    return hline
   
plt.sca(axHist)
yLims = [0, 15]
binEdges = np.linspace(-1, 1, 32) #32
for indType, cellType in enumerate(cellTypes):
    colorThisType = figparams.colors[cellType]
    hline = line_histogram(toneRespIndexEachType[indType], binEdges, lw=2,
                           percent=True, color=colorThisType)
    #hline = line_histogram(np.abs(toneRespIndexEachType[indType]), binEdges, lw=2,
    #                       percent=True, color=colorThisType)
    medianPos = np.median(toneRespIndexEachTypePos[indType])
    medianNeg = np.median(toneRespIndexEachTypeNeg[indType])
    plt.plot(medianPos, 0.95*yLims[-1], 'v', color=colorThisType, mew=2, mfc='none')
    plt.plot(medianNeg, 0.95*yLims[-1], 'v', color=colorThisType, mew=2, mfc='none')
    
plt.ylim(yLims)
#plt.ylim([0,60])
extraplots.boxoff(axHist)
plt.yticks(np.arange(0,20,5))
plt.xlabel('Tone response index')
plt.ylabel('% cells')
plt.axvline(0, color='0.9')
#plt.legend(['D1','non-D1'], loc='upper left')
plt.text(-1, 11, 'D1 cells', ha='left', fontweight='bold', fontsize=fontSizeLabels+0,
         color=figparams.colors['D1'])
plt.text(-1, 9, 'non-D1 cells', ha='left', fontweight='bold', fontsize=fontSizeLabels+0,
         color=figparams.colors['ND1'])

plt.show()


uval, pValPos = stats.mannwhitneyu(toneRespIndexD1pos, toneRespIndexND1pos, alternative='two-sided')
print(f'Tone index (pos): D1: {toneRespIndexD1pos.median():0.4f} vs ' +
      f'ND1: {toneRespIndexND1pos.median():0.4f}   ' +
      f'p = {pValPos:0.6f} (U={uval})')
uval, pValNeg = stats.mannwhitneyu(toneRespIndexD1neg, toneRespIndexND1neg, alternative='two-sided')
print(f'Tone index (neg): D1: {toneRespIndexD1neg.median():0.4f} vs ' +
      f'ND1: {toneRespIndexND1neg.median():0.4f}   ' +
      f'p = {pValNeg:0.6f} (U={uval})')
uval, pValNeg = stats.mannwhitneyu(np.abs(toneRespIndexD1[~np.isnan(toneRespIndexD1)]),
                                   np.abs(toneRespIndexND1[~np.isnan(toneRespIndexND1)]),
                                   alternative='two-sided')
print(f'Tone index (abs): ' +
      f'D1 ({len(toneRespIndexD1)}): {np.abs(toneRespIndexD1).median():0.4f} vs ' +
      f'ND1 ({len(toneRespIndexND1)}): {np.abs(toneRespIndexND1).median():0.4f}   ' +
      f'p = {pValNeg:0.6f} (U={uval})')
print()


print(f'Tone responsive: D1: {np.sum(indD1 & toneResponsive)} ' +
      f'ND1: {np.sum(indND1 & toneResponsive)} ')
print()


# -------------------- Bandwidth ------------------------
goodFit = celldb.toneGaussianRsquare > 0.01 #0.01
#goodFit = (celldb.toneGaussianRsquare > 0.01) & (celldb.toneGaussianX0<15) & (celldb.toneGaussianX0>11)
fullWidthHalfMax = 2.355 * celldb.toneGaussianSigma
bwD1 = fullWidthHalfMax[indD1 & toneResponsive & goodFit]
bwND1 = fullWidthHalfMax[indND1 & toneResponsive & goodFit]
nD1bw = len(bwD1)
nND1bw = len(bwND1)
medianD1bw = np.median(bwD1)
medianND1bw = np.median(bwND1)

plt.sca(axBW)
markerSize = 3
markerWidth = 0.5
lineWidth = 3
jitter = 0.2
xLine = 0.3*np.array([-1,1])
np.random.seed(0)
xpos = 1 + jitter*(2*np.random.rand(nD1bw)-1)
plt.plot(xpos, bwD1, 'o', mfc='none', mec=figparams.colors['D1'], ms=markerSize, mew=markerWidth)
plt.plot(1+xLine, [medianD1bw, medianD1bw], lw=lineWidth, color=figparams.colors['D1'])
xpos = 2 + jitter*(2*np.random.rand(nND1bw)-1)
plt.plot(xpos, bwND1, 'o', mfc='none', mec=figparams.colors['ND1'], ms=markerSize, mew=markerWidth)
plt.plot(2+xLine, [medianND1bw, medianND1bw], lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axBW.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axBW)
plt.ylabel('Tuning bandwidth (oct)', fontsize=fontSizeLabels)

uval, pValBW = stats.mannwhitneyu(bwD1, bwND1, alternative='two-sided')
print(f'Tuning BW: D1 ({nD1bw}): {medianD1bw:0.4f} vs ' +
      f'ND1 ({nND1bw}): {medianND1bw:0.4f}   ' +
      f'p = {pValBW:0.6f} (U={uval})\n')


# -------------------- Onset vs Sustain ------------------------
baselineFiring = celldb['toneFiringRateBaseline']
onsetResp = celldb['toneFiringRateBestOnset']
sustainResp = celldb['toneFiringRateBestSustain']
onsetToSustainIndex = (onsetResp - sustainResp) / (onsetResp + sustainResp)

toneResponsive = toneResponsive & goodFit; print('WARNING! only good fit neurons used here!')
#toneResponsive = toneResponsive & ((celldb['toneFiringRateBest']-celldb['toneFiringRateBaseline'])>0)

# NOTE: I need to remove nans because some neurons may be have FR=0 for tones but not AM
osiD1 = onsetToSustainIndex[indD1 & toneResponsive & ~np.isnan(onsetToSustainIndex)]
osiND1 = onsetToSustainIndex[indND1 & toneResponsive & ~np.isnan(onsetToSustainIndex)]
nD1osi = len(osiD1)
nND1osi = len(osiND1)
medianD1osi = np.median(osiD1)
medianND1osi = np.median(osiND1)

'''
########3# TEST #########
np.sum(np.isnan(osiD1))
np.sum(np.isnan(onsetToSustainIndex))
np.sum(np.isnan(onsetToSustainIndex[toneResponsive]))
cc = celldb[toneResponsive]
inds = np.flatnonzero(np.isnan(onsetToSustainIndex[toneResponsive]))
print(inds) #  [147, 176, 204]
cc.iloc[147]
# It looks like some with low FR are passing criteria in studyutils.select_cells()
thFR = 1.0
highFR = ( (celldb.toneFiringRateBaseline > thFR) | (celldb.toneFiringRateBest > thFR) |
           (celldb.amFiringRateBaseline > thFR) | (celldb.amFiringRateBestOnset > thFR) |
           (celldb.amFiringRateBestSustain > thFR))
celldb[~highFR & indD1]
celldb[~highFR & toneResponsive]
# look at 3015 it seems to have all FR low but be tone responsive
#######################3
'''

plt.sca(axOSI)
#markerSize = 3
#lineWidth = 3
#jitter = 0.2
#xLine = 0.3*np.array([-1,1])
#np.random.seed(0)
xpos = 1 + jitter*(2*np.random.rand(nD1osi)-1)
plt.plot(xpos, osiD1, 'o', mfc='none', mec=figparams.colors['D1'], ms=markerSize, mew=markerWidth)
plt.plot(1+xLine, [medianD1osi, medianD1osi], lw=lineWidth, color=figparams.colors['D1'])
xpos = 2 + jitter*(2*np.random.rand(nND1osi)-1)
plt.plot(xpos, osiND1, 'o', mfc='none', mec=figparams.colors['ND1'], ms=markerSize, mew=markerWidth)
plt.plot(2+xLine, [medianND1osi, medianND1osi], lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axOSI.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axOSI)
plt.ylabel('Onset-to-sustain index', fontsize=fontSizeLabels)

uval, pValOSI = stats.mannwhitneyu(osiD1, osiND1, alternative='two-sided')
print(f'Onset-to-sustain index: D1 ({nD1osi}): {medianD1osi:0.4f} vs ' +
      f'ND1 ({nND1osi}): {medianND1osi:0.4f}   ' +
      f'p = {pValOSI:0.6f} (U={uval})\n')




plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
