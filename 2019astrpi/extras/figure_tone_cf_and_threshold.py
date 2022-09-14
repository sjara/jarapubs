"""
Figure about responses to pure tones.
"""

import sys
sys.path.append('..')
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

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_tone_cf_and_threshold' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [4.6, 3.8] # In inches
figSize = [7, 3.6] # In inches


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

#panelLabels = ['A', ['B','D'], ['C','E']]
#labelPosX = [0.04, 0.3, 0.68]   # Horiz position for panel labels
labelPosXtop = [0.01, 0.30, 0.525, 0.755]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.55, 0.745]   # Horiz position for panel labels
labelPosY = [0.93, 0.5]    # Vert position for panel labels

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220212.h5')
#dbPath = os.path.join(figuresDataDir, 'astrpi_am_tuning.h5')
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)

pp = lambda x: f"['{x.subject}', '{x.date}', {x.pdepth}, {x.egroup}, {x.cluster}]"
# pp(celldb.loc[])

# D1: 1967, 1971
selectedD1 = ['d1pi036', '2019-05-29', 2800, 2, 4] # cell 710
selectedND1 = ['d1pi041', '2019-08-25', 3600, 7, 2] # cell 1454
selectedD1sup = ['d1pi046', '2020-02-18', 3200, 8, 3] # cell 2999
selectedND1sup = ['d1pi041', '2019-08-25', 3400, 6, 2] # cell 1381

#selectedND1sup = ['d1pi041', '2019-08-25', 3600, 5, 4] # cell 1446
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 2, 3]  # cell 1971
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 1, 2]  # cell 1967
#selectedND1 = ['d1pi041', '2019-08-25', 3500, 8, 2] # cell 1426

cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
selectedCells = [selectedD1, selectedND1, selectedD1sup, selectedND1sup]
examplesTuning = ['D1', 'ND1', 'D1', 'ND1']
cMaps = ['Blues', 'Reds', 'Blues', 'Reds']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 1)
gsMain.update(left=0.09, right=0.95, top=1.0, bottom=0.12, hspace=0.3)

gsTC = gsMain[0,0].subgridspec(1, 4, wspace=0.3)
axTCD1 = plt.subplot(gsTC[0, 0])
axTCD1.annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCND1 = plt.subplot(gsTC[0, 1])
axTCND1.annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
axTCSupD1 = plt.subplot(gsTC[0, 2])
axTCSupD1.annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axTCSupND1 = plt.subplot(gsTC[0, 3])
axTCSupND1.annotate('D', xy=(labelPosXtop[3],labelPosY[0]), xycoords='figure fraction',
                    fontsize=fontSizePanel, fontweight='bold')
axTC = [axTCD1, axTCND1, axTCSupD1, axTCSupND1]

gsBottom = gsMain[1,0].subgridspec(1, 3, wspace=0.6, width_ratios=[0.6, 0.15, 0.15])
axHist = plt.subplot(gsBottom[0, 0])
axHist.annotate('E', xy=(labelPosXbot[0],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axThresh = plt.subplot(gsBottom[0, 1])
axThresh.annotate('F', xy=(labelPosXbot[1],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axCF = plt.subplot(gsBottom[0, 2])
axCF.annotate('G', xy=(labelPosXbot[2],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')


if 0:
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
        if indType in [2,3]:
            responsePeriod = [0, 0.2];
            print('WARNING! Using 0.2s period for this example.\n')
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
        if indType==2: # This cell used to be #2 (before including D1sup)
            maxFR = 10; print('WARNING! fixing max color for third example.\n')
        if indType==3: # This cell used to be #2 (before including D1sup)
            maxFR = 24; print('WARNING! fixing max color for fourth example.\n')
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
        else:
            axTC[indType].set_yticklabels([])
        extraplots.set_ticks_fontsize(axTC[indType], fontSizeTicks)

        cbar = plt.colorbar(tuningPlot, ax=axTC[indType], format='%d', shrink=0.52,
                            pad=0.06, aspect=10)
        cbar.ax.set_ylabel('spk/s', fontsize=fontSizeLabels, labelpad=-2, rotation=-90)
        extraplots.set_ticks_fontsize(cbar.ax, fontSizeTicks)
        cbar.set_ticks([0, maxFR])
        cbar.ax.axhline(dbRow.toneFiringRateBaseline, color=[0,1,0])
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
    baselineD1 = celldb.toneFiringRateBaseline[indD1 & cellsWithTone]
    baselineND1 = celldb.toneFiringRateBaseline[indND1 & cellsWithTone]

toneRespIndexD1pos = toneRespIndexD1[toneRespIndexD1>0]
toneRespIndexND1pos = toneRespIndexND1[toneRespIndexND1>0]
toneRespIndexD1neg = toneRespIndexD1[toneRespIndexD1<0]
toneRespIndexND1neg = toneRespIndexND1[toneRespIndexND1<0]
toneRespIndexEachTypePos = [toneRespIndexD1pos, toneRespIndexND1pos]
toneRespIndexEachTypeNeg = [toneRespIndexD1neg, toneRespIndexND1neg]
baselineD1pos = baselineD1[toneRespIndexD1>0]
baselineD1neg = baselineD1[toneRespIndexD1<0]
baselineND1pos = baselineND1[toneRespIndexND1>0]
baselineND1neg = baselineND1[toneRespIndexND1<0]


def line_histogram(data, edges, percent=False, **kwargs):
    hh, ee = np.histogram(data, edges)
    if percent:
        hh = 100 * hh/len(data[~np.isnan(data)])
    #print(np.sum(hh))
    hline = plt.step(np.r_[ee, ee[-1]], np.r_[0,hh,0],'-', where='pre', **kwargs)
    return hline

'''
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
'''

nPosD1 = len(toneRespIndexD1pos)
nNegD1 = len(toneRespIndexD1neg)
nTotalD1 = nPosD1+nNegD1
nPosND1 = len(toneRespIndexND1pos)
nNegND1 = len(toneRespIndexND1neg)
nTotalND1 = nPosND1+nNegND1
print(f'Neg responses (D1): {nNegD1}/{nTotalD1} = {nNegD1/nTotalD1:0.01%}')
print(f'Neg responses (ND1): {nNegND1}/{nTotalND1} = {nNegND1/nTotalND1:0.01%}')

uval, pValBaselineD1 = stats.mannwhitneyu(baselineD1pos, baselineD1neg, alternative='two-sided')
print(f'Baseline firing D1: pos: {baselineD1pos.median():0.4f} vs ' +
      f'neg: {baselineD1neg.median():0.4f}   ' +
      f'p = {pValBaselineD1:0.6f} (U={uval})')

uval, pValBaselineND1 = stats.mannwhitneyu(baselineND1pos, baselineND1neg, alternative='two-sided')
print(f'Baseline firing ND1: pos: {baselineND1pos.median():0.4f} vs ' +
      f'neg: {baselineND1neg.median():0.4f}   ' +
      f'p = {pValBaselineND1:0.6f} (U={uval})')
uval, pValBaselinePos = stats.mannwhitneyu(baselineD1pos, baselineND1pos, alternative='two-sided')
print(f'Baseline firing (pos): D1: {baselineD1pos.median():0.4f} vs ' +
      f'neg: {baselineND1pos.median():0.4f}   ' +
      f'p = {pValBaselinePos:0.6f} (U={uval})')
uval, pValBaselineNeg = stats.mannwhitneyu(baselineD1neg, baselineND1neg, alternative='two-sided')
print(f'Baseline firing (neg): D1: {baselineD1neg.median():0.4f} vs ' +
      f'neg: {baselineND1neg.median():0.4f}   ' +
      f'p = {pValBaselineNeg:0.6f} (U={uval})')
uval, pValBaseline = stats.mannwhitneyu(baselineD1, baselineND1, alternative='two-sided')
print(f'Baseline firing: D1: {baselineD1.median():0.4f} vs ' +
      f'neg: {baselineND1.median():0.4f}   ' +
      f'p = {pValBaseline:0.6f} (U={uval})')
print()

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


# -------------------- Threshold ------------------------
goodFit = celldb.toneGaussianRsquare > 0.01 #0.01
#goodFit = (celldb.toneGaussianRsquare > 0.01) & (celldb.toneGaussianX0<15) & (celldb.toneGaussianX0>11)
#fullWidthHalfMax = 2.355 * celldb.toneGaussianSigma

intensityThreshold = celldb.toneIntensityThreshold
validThreshold = intensityThreshold<np.inf  #65

threshD1 = intensityThreshold[indD1 & toneResponsive & goodFit & validThreshold]
threshND1 = intensityThreshold[indND1 & toneResponsive & goodFit & validThreshold]
#threshD1 = intensityThreshold[indD1 & toneResponsive & validThreshold]
#threshND1 = intensityThreshold[indND1 & toneResponsive & validThreshold]
nD1thresh = len(threshD1)
nND1thresh = len(threshND1)
medianD1thresh = np.median(threshD1)
medianND1thresh = np.median(threshND1)
###medianD1thresh = np.mean(threshD1)
###medianND1thresh = np.mean(threshND1)

plt.sca(axThresh)
markerSize = 3
markerWidth = 0.5
lineWidth = 3
jitter = 0.2
xLine = 0.3*np.array([-1,1])
np.random.seed(0)
xpos = 1 + jitter*(2*np.random.rand(nD1thresh)-1)
plt.plot(xpos, threshD1, 'o', mfc='none', mec=figparams.colors['D1'], ms=markerSize, mew=markerWidth)
#plt.plot(1+xLine, [medianD1thresh, medianD1thresh], lw=lineWidth, color=figparams.colors['D1'])
xpos = 2 + jitter*(2*np.random.rand(nND1thresh)-1)
plt.plot(xpos, threshND1, 'o', mfc='none', mec=figparams.colors['ND1'], ms=markerSize, mew=markerWidth)
#plt.plot(2+xLine, [medianND1thresh, medianND1thresh], lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axThresh.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axThresh)
plt.ylabel('Intensity threshold (dB)', fontsize=fontSizeLabels)

uval, pValThresh = stats.mannwhitneyu(threshD1, threshND1, alternative='two-sided')
print(f'Tuning Thresh: D1 ({nD1thresh}): {medianD1thresh:0.4f} vs ' +
      f'ND1 ({nND1thresh}): {medianND1thresh:0.4f}   ' +
      f'p = {pValThresh:0.6f} (U={uval})\n')


# -------------------- CF ------------------------
goodFit = celldb.toneGaussianRsquare > 0.01 #0.01
#goodFit = (celldb.toneGaussianRsquare > 0.01) & (celldb.toneGaussianX0<15) & (celldb.toneGaussianX0>11)
#fullWidthHalfMax = 2.355 * celldb.toneGaussianSigma

charactFreq = celldb.toneCharactFreq
###charactFreq = 2**celldb.toneGaussianX0; print('WARNING! Using center of Gaussian, not CF')

cfD1 = charactFreq[indD1 & toneResponsive & goodFit]
cfND1 = charactFreq[indND1 & toneResponsive & goodFit]
nD1cf = len(cfD1)
nND1cf = len(cfND1)
medianD1cf = np.median(cfD1)
medianND1cf = np.median(cfND1)

plt.sca(axCF)
markerSize = 3
markerWidth = 0.5
lineWidth = 3
jitter = 0.2
xLine = 0.3*np.array([-1,1])
np.random.seed(0)
xpos = 1 + jitter*(2*np.random.rand(nD1cf)-1)
plt.plot(xpos, cfD1, 'o', mfc='none', mec=figparams.colors['D1'], ms=markerSize, mew=markerWidth)
plt.plot(1+xLine, [medianD1cf, medianD1cf], lw=lineWidth, color=figparams.colors['D1'])
xpos = 2 + jitter*(2*np.random.rand(nND1cf)-1)
plt.plot(xpos, cfND1, 'o', mfc='none', mec=figparams.colors['ND1'], ms=markerSize, mew=markerWidth)
plt.plot(2+xLine, [medianND1cf, medianND1cf], lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axCF.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axCF)
plt.ylabel('CF (Hz)', fontsize=fontSizeLabels)

uval, pValCF = stats.mannwhitneyu(cfD1, cfND1, alternative='two-sided')
print(f'Tuning CF: D1 ({nD1cf}): {medianD1cf:0.4f} vs ' +
      f'ND1 ({nND1cf}): {medianND1cf:0.4f}   ' +
      f'p = {pValCF:0.6f} (U={uval})\n')


# -------------------- Threshold vs CF ------------------------
plt.sca(axHist)
from matplotlib.ticker import ScalarFormatter
nCells = len(charactFreq)
xjitter = 0.05
yjitter = 1
xpos = 1 + xjitter*(2*np.random.rand(nCells)-1)
ypos = yjitter*(2*np.random.rand(nCells)-1)
plt.plot(charactFreq*xpos, intensityThreshold+ypos,'o', mfc='none')
axHist.set_xscale('log')
axHist.xaxis.set_major_formatter(ScalarFormatter())
plt.xlabel('Frequency (Hz)')
plt.ylabel('Intensity threshold (dB)')

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
