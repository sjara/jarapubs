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

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'fig3_tone_responses' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [4.6, 3.8] # In inches
figSize = [7, 4.5] # In inches

PANELS = [1, 1, 1,1,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

#panelLabels = ['A', ['B','D'], ['C','E']]
#labelPosX = [0.04, 0.3, 0.68]   # Horiz position for panel labels
labelPosXtop = [0.01, 0.30, 0.535, 0.765]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.43, 0.63, 0.8]   # Horiz position for panel labels
labelPosY = [0.95, 0.72, 0.38]    # Vert position for panel labels

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220212.h5')
#dbPath = os.path.join(figuresDataDir, 'astrpi_am_tuning.h5')
dbPath = os.path.join(figuresDataDir, 'astrpi_freq_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)

pp = lambda x: f"['{x.subject}', '{x.date}', {x.pdepth}, {x.egroup}, {x.cluster}]"
# pp(celldb.loc[])

# D1: 1967, 1971
selectedD1sym = ['d1pi036', '2019-05-29', 2800, 2, 4] # cell 710
selectedD1asym = ['d1pi042', '2019-09-11', 3300, 2, 3] # cell 1971
selectedD1sup = ['d1pi036', '2019-05-29', 2800, 4, 4] # cell 717
selectedD1ht = ['d1pi047', '2020-02-24', 3600, 6, 2] # cell 3927

selectedND1sym = ['d1pi041', '2019-08-25', 3600, 7, 2] # cell 1454
selectedND1asym = ['d1pi041', '2019-08-25', 3400, 7, 6] # cell 1391
selectedND1sup = ['d1pi041', '2019-08-25', 3400, 6, 2] # cell 1381
selectedND1ht = ['d1pi049', '2020-03-11', 3400, 8, 5] # cell 4336

'''
# Testing high-threshold cells
selectedD1sym = ['d1pi032', '2019-02-22', 2900, 8, 2]
selectedD1asym =['d1pi041', '2019-08-25', 2700, 8, 5]
selectedD1ht = ['d1pi049', '2020-03-07', 3400, 8, 4]
selectedD1sup = ['d1pi047', '2020-02-24', 3600, 6, 2]
selectedND1sym = ['d1pi049', '2020-03-11', 3400, 8, 5]
selectedND1asym =['d1pi041', '2019-08-25', 3100, 7, 3]
selectedND1ht = ['d1pi041', '2019-08-27', 3400, 4, 3]
selectedND1sup = ['d1pi047', '2020-02-24', 3400, 4, 5]
'''

#selectedND1sup = ['d1pi041', '2019-08-25', 3600, 5, 4] # cell 1446
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 2, 3]  # cell 1971
#selectedD1 = ['d1pi042', '2019-09-11', 3300, 1, 2]  # cell 1967
#selectedND1 = ['d1pi041', '2019-08-25', 3500, 8, 2] # cell 1426
#selectedND1asym = ['d1pi046', '2020-02-16', 3600, 6, 5] # cell 2869
#selectedD1sup = ['d1pi046', '2020-02-18', 3200, 8, 3] # cell 2999

cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
selectedCells = [selectedD1sym, selectedD1asym, selectedD1sup, selectedD1ht,
                 selectedND1sym, selectedND1asym, selectedND1sup, selectedND1ht]
examplesTuning = ['D1', 'D1', 'D1', 'D1', 'ND1', 'ND1', 'ND1', 'ND1']
cMaps = ['Blues', 'Blues', 'Blues', 'Blues', 'Reds', 'Reds', 'Reds', 'Reds']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 1, height_ratios=[0.6, 0.4])
gsMain.update(left=0.09, right=0.97, top=0.97, bottom=0.1, hspace=0.4)

gsTCD1 = gsMain[0,0].subgridspec(2, 4, wspace=0.3)
axTCD1sym = plt.subplot(gsTCD1[0, 0])
axTCD1sym.annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCD1asym = plt.subplot(gsTCD1[0, 1])
axTCD1sym.annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCD1ht = plt.subplot(gsTCD1[0, 2])
axTCD1ht.annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCD1sup = plt.subplot(gsTCD1[0, 3])
axTCD1sup.annotate('D', xy=(labelPosXtop[3],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
#gsTCND1 = gsMain[1,0].subgridspec(1, 4, wspace=0.3)
axTCND1sym = plt.subplot(gsTCD1[1, 0])
axTCND1sym.annotate('E', xy=(labelPosXtop[0],labelPosY[1]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
axTCND1asym = plt.subplot(gsTCD1[1, 1])
axTCND1sym.annotate('F', xy=(labelPosXtop[1],labelPosY[1]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')
axTCND1ht = plt.subplot(gsTCD1[1, 2])
axTCND1ht.annotate('G', xy=(labelPosXtop[2],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axTCND1sup = plt.subplot(gsTCD1[1, 3])
axTCND1sup.annotate('H', xy=(labelPosXtop[3],labelPosY[1]), xycoords='figure fraction',
                    fontsize=fontSizePanel, fontweight='bold')
axTC = [axTCD1sym, axTCD1asym, axTCD1ht, axTCD1sup,
        axTCND1sym, axTCND1asym, axTCND1ht, axTCND1sup]

gsBottom = gsMain[1,0].subgridspec(1, 4, wspace=0.7, width_ratios=[0.6, 0.15, 0.15, 0.15])
axHist = plt.subplot(gsBottom[0, 0])
axHist.annotate('I', xy=(labelPosXbot[0],labelPosY[2]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axThresh = plt.subplot(gsBottom[0, 1])
axThresh.annotate('J', xy=(labelPosXbot[1],labelPosY[2]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axBW = plt.subplot(gsBottom[0, 2])
axBW.annotate('K', xy=(labelPosXbot[2],labelPosY[2]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axFRAtypes = plt.subplot(gsBottom[0, 3])
axFRAtypes.annotate('L', xy=(labelPosXbot[3],labelPosY[2]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
#axOSI = plt.subplot(gsBottom[0, 3])
#axOSI.annotate('H', xy=(labelPosXbot[2],labelPosY[1]), xycoords='figure fraction',
#                fontsize=fontSizePanel, fontweight='bold')


if PANELS[0]:
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
        '''
        if indType in [2,5]:
            responsePeriod = [0, 0.2];
            print('WARNING! Using 0.2s period for this example.\n')
        '''
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
        if indRow==2999: # This cell used to be #2 (before including D1sup)
            maxFR = 8; print('WARNING! fixing max color for suppression example.\n')
        #if indRow==1381: # This cell used to be #2 (before including D1sup)
        #    maxFR = 24; print('WARNING! fixing max color for suppresion example.\n')
        if indRow==4336000:
            maxFR = 20; print('WARNING! fixing max color for suppression example.\n')
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
        if indType > 3:
            axTC[indType].set_xticklabels(freqLabels)
            plt.xlabel('Frequency (kHz)', fontsize=fontSizeLabels)
        else:
            axTC[indType].set_xticklabels([])
        if indType in [0, 4]:
            plt.ylabel('Intensity\n(dB SPL)', fontsize=fontSizeLabels)
        else:
            axTC[indType].set_yticklabels([])
        extraplots.set_ticks_fontsize(axTC[indType], fontSizeTicks)

        # shrink = 0.52
        cbar = plt.colorbar(tuningPlot, ax=axTC[indType], format='%d', shrink=0.84,
                            pad=0.06, aspect=10)
        cbar.ax.set_ylabel('spk/s', fontsize=fontSizeLabels, labelpad=-2, rotation=-90)
        extraplots.set_ticks_fontsize(cbar.ax, fontSizeTicks)
        cbar.set_ticks([0, maxFR])
        cbar.ax.axhline(dbRow.toneFiringRateBaseline, color=[0,0,0])
        plt.show()
        #sys.exit()


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


if PANELS[1]:
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

nPosD1 = len(toneRespIndexD1pos)
nNegD1 = len(toneRespIndexD1neg)
nTotalD1 = nPosD1+nNegD1
nPosND1 = len(toneRespIndexND1pos)
nNegND1 = len(toneRespIndexND1neg)
nTotalND1 = nPosND1+nNegND1

negRespTable = [ [nNegD1, nNegND1],
                         [nTotalD1-nNegD1, nTotalND1-nNegND1] ]
oddsratio, pValNegResp = stats.fisher_exact(negRespTable, alternative='two-sided')
print(f'Neg responses (D1): {nNegD1}/{nTotalD1} = {nNegD1/nTotalD1:0.01%}')
print(f'Neg responses (ND1): {nNegND1}/{nTotalND1} = {nNegND1/nTotalND1:0.01%}')
print(f' p = {pValNegResp:0.6f}  (oddratio={oddsratio})')
print('')

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
print(f'Baseline firing: D1: {baselineD1.median():0.4f} ({len(baselineD1)}) vs ' +
      f'ND1: {baselineND1.median():0.4f} ({len(baselineND1)})   ' +
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
goodFit = celldb.toneGaussianRsquare > 0.01

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

uval, pValThresh = stats.mannwhitneyu(threshD1, threshND1, alternative='two-sided')
print(f'Tuning Thresh: D1 ({nD1thresh}): {medianD1thresh:0.4f} vs ' +
      f'ND1 ({nND1thresh}): {medianND1thresh:0.4f}   ' +
      f'p = {pValThresh:0.6f} (U={uval})\n')

# --- Characteristic frequency ---
charactFreq = celldb.toneCharactFreq
validCF = ~np.isnan(charactFreq)
cfD1 = charactFreq[indD1 & toneResponsive & goodFit & validCF]
cfND1 = charactFreq[indND1 & toneResponsive & goodFit & validCF]
#cfD1 = charactFreq[indD1 & toneResponsive  & validCF]
#cfND1 = charactFreq[indND1 & toneResponsive  & validCF]
nD1cf = len(cfD1)
nND1cf = len(cfND1)
medianD1cf = np.median(cfD1)
medianND1cf = np.median(cfND1)
uval, pValCF = stats.mannwhitneyu(cfD1, cfND1, alternative='two-sided')
print(f'Tuning Charact Freq: D1 ({nD1cf}): {medianD1cf:0.4f} vs ' +
      f'ND1 ({nND1cf}): {medianND1cf:0.4f}   ' +
      f'p = {pValCF:0.6f} (U={uval})')
print(f'    CF: D1  range = [{cfD1.min()}, {cfD1.max()}]')
print(f'    CF: ND1 range = [{cfND1.min()}, {cfND1.max()}]')
print('')

# --- Best frequency ---
bestFreqInd = celldb.toneIndexBest
validBF = bestFreqInd>0
bfD1 = possibleFreq[bestFreqInd[indD1 & toneResponsive & goodFit]]
bfND1 = possibleFreq[bestFreqInd[indND1 & toneResponsive & goodFit]]
nD1bf = len(bfD1)
nND1bf = len(bfND1)
medianD1bf = np.median(bfD1)
medianND1bf = np.median(bfND1)
uval, pValBF = stats.mannwhitneyu(bfD1, bfND1, alternative='two-sided')
print(f'Tuning Best Freq: D1 ({nD1bf}): {medianD1bf:0.4f} vs ' +
      f'ND1 ({nND1bf}): {medianND1bf:0.4f}   ' +
      f'p = {pValBF:0.6f} (U={uval})')
print(f'    BF: D1  range = [{bfD1.min()}, {bfD1.max()}]')
print(f'    BF: ND1 range = [{bfND1.min()}, {bfND1.max()}]')
print('')


plt.sca(axThresh)
if 1:
    boxWidth = 0.5
    lineWidth = 1.2
    pboxD1 = plt.boxplot(threshD1, positions=[1], widths=[boxWidth])
    pboxND1 = plt.boxplot(threshND1, positions=[2], widths=[boxWidth])
    plt.setp(pboxD1['boxes'][0], color=figparams.colors['D1'], lw=lineWidth)
    plt.setp(pboxD1['medians'][0], color=figparams.colors['D1'], lw=lineWidth)
    plt.setp(pboxD1['whiskers'], color=figparams.colors['D1'])
    plt.setp(pboxD1['caps'], color=figparams.colors['D1'])
    plt.setp(pboxND1['boxes'][0], color=figparams.colors['ND1'], lw=lineWidth)
    plt.setp(pboxND1['medians'][0], color=figparams.colors['ND1'], lw=lineWidth)
    plt.setp(pboxND1['whiskers'], color=figparams.colors['ND1'])
    plt.setp(pboxND1['caps'], color=figparams.colors['ND1'])
    plt.ylim([10,75])
    if pValThresh<0.05:
        hs, hl = extraplots.significance_stars([1, 2], 75, yLength=0.9, gapFactor=0.3, starSize=6)
        plt.setp(hl, lw=0.75)

    '''
    vpD1 = plt.violinplot(threshD1, positions=[1], showextrema=False, showmedians=True)
    vpD1['bodies'][0].set_facecolor(figparams.colors['D1'])
    vpD1['cmedians'].set_color(figparams.colors['D1'])
    vpD1 = plt.violinplot(threshND1, positions=[2], showextrema=False, showmedians=True)
    vpD1['bodies'][0].set_facecolor(figparams.colors['ND1'])
    vpD1['cmedians'].set_color(figparams.colors['ND1'])
    '''
else:
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
jitter = 0.2
if 1:
    boxWidth = 0.5
    lineWidth = 1.2
    pboxD1 = plt.boxplot(bwD1, positions=[1], widths=[boxWidth])
    pboxND1 = plt.boxplot(bwND1, positions=[2], widths=[boxWidth])
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
    plt.ylim([0, 2.5])
    '''
    vpD1 = plt.violinplot(bwD1, positions=[1], showextrema=False, showmedians=True)
    vpD1['bodies'][0].set_facecolor(figparams.colors['D1'])
    vpD1['cmedians'].set_color(figparams.colors['D1'])
    vpD1 = plt.violinplot(bwND1, positions=[2], showextrema=False, showmedians=True)
    vpD1['bodies'][0].set_facecolor(figparams.colors['ND1'])
    vpD1['cmedians'].set_color(figparams.colors['ND1'])
    '''
else:
    plt.sca(axBW)
    markerSize = 3
    markerWidth = 0.5
    lineWidth = 3
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
      f'p = {pValBW:0.6f} (U={uval})')

# -- Additional bandwidth analysis (BW10 and BW40) --
haveBW10 = ~np.isnan(celldb.toneBW10)
haveBW40 = ~np.isnan(celldb.toneBW40)
haveCF = ~np.isnan(celldb.toneCharactFreq)
validFRA = toneResponsive & haveCF
BW10D1 = celldb.toneBW10[indD1 & validFRA & haveBW10]
BW10ND1 = celldb.toneBW10[indND1 & validFRA & haveBW10]
uval, pValBW10 = stats.mannwhitneyu(BW10D1, BW10ND1, alternative='two-sided')
print(f'Tuning BW10: D1 ({len(BW10D1)}): {np.median(BW10D1):0.4f} vs ' +
      f'ND1 ({len(BW10ND1)}): {np.median(BW10ND1):0.4f}   ' +
      f'p = {pValBW10:0.6f} (U={uval})')
BW40D1 = celldb.toneBW40[indD1 & validFRA & haveBW40]
BW40ND1 = celldb.toneBW40[indND1 & validFRA & haveBW40]
uval, pValBW40 = stats.mannwhitneyu(BW40D1, BW40ND1, alternative='two-sided')
print(f'Tuning BW40: D1 ({len(BW40D1)}): {np.median(BW40D1):0.4f} vs ' +
      f'ND1 ({len(BW40ND1)}): {np.median(BW40ND1):0.4f}   ' +
      f'p = {pValBW40:0.6f} (U={uval})')
print('')


# -------------------- Types of FRA ------------------------
haveCF = ~np.isnan(celldb.toneCharactFreq)
lowThreshold = intensityThreshold<60
negResponse = celldb.toneGaussianA<0
validFRA = toneResponsive & haveCF & lowThreshold & ~negResponse
validFRAneg = toneResponsive & haveCF & lowThreshold & negResponse
lowSlopeD1 = celldb.toneFRAlowSlope[indD1 & validFRA]
highSlopeD1 = celldb.toneFRAhighSlope[indD1 & validFRA]
lowSlopeND1 = celldb.toneFRAlowSlope[indND1 & validFRA]
highSlopeND1 = celldb.toneFRAhighSlope[indND1 & validFRA]
highSlopeThreshold = 100
lowSlopeThreshold = 100

symmetricTuningD1 = (highSlopeD1<highSlopeThreshold)
asymmetricTuningD1 = (highSlopeD1>=highSlopeThreshold)
symmetricTuningND1 = (highSlopeND1<highSlopeThreshold)
asymmetricTuningND1 = (highSlopeND1>=highSlopeThreshold)
'''
symmetricTuningD1 = (highSlopeD1<highSlopeThreshold) & (-lowSlopeD1<lowSlopeThreshold)
asymmetricTuningD1 = ~symmetricTuningD1
symmetricTuningND1 = (highSlopeND1<highSlopeThreshold) & (-lowSlopeND1<lowSlopeThreshold)
asymmetricTuningND1 = ~symmetricTuningND1
'''
'''
symmetricTuningD1 = (-lowSlopeD1>lowSlopeThreshold)
asymmetricTuningD1 = ~symmetricTuningD1
symmetricTuningND1 = (-lowSlopeND1>lowSlopeThreshold)
asymmetricTuningND1 = ~symmetricTuningND1
'''
'''
symmetricTuningD1 = (highSlopeD1<highSlopeThreshold) & (-lowSlopeD1<highSlopeThreshold)
asymmetricTuningD1 = (highSlopeD1>=highSlopeThreshold)
symmetricTuningND1 = (highSlopeND1<highSlopeThreshold) & (-lowSlopeND1<highSlopeThreshold)
asymmetricTuningND1 = (highSlopeND1>=highSlopeThreshold)
'''
nD1neg = np.sum(validFRAneg & indD1)
nND1neg = np.sum(validFRAneg & indND1)
nD1 = len(symmetricTuningD1) + nD1neg
nND1 = len(symmetricTuningND1) + nND1neg
nSymmetricTuningD1 = np.sum(symmetricTuningD1)
nAsymmetricTuningD1 = np.sum(asymmetricTuningD1)
nSymmetricTuningND1 = np.sum(symmetricTuningND1)
nAsymmetricTuningND1 = np.sum(asymmetricTuningND1)

print('---------- Types of FRA ----------')
cTable1 = [[sum(symmetricTuningD1), sum(symmetricTuningND1)],
          [sum(asymmetricTuningD1), sum(asymmetricTuningND1)]]
cTable2 = [[nAsymmetricTuningD1, nAsymmetricTuningND1],
           [nD1-nAsymmetricTuningD1, nD1-nAsymmetricTuningND1]]
oddsr, pvalFRAtypes = stats.fisher_exact(cTable2)  
print(f'D1: sym = {sum(symmetricTuningD1)/nD1:0.2%} ({sum(symmetricTuningD1)}/{nD1}), ' + \
      f'asym = {sum(asymmetricTuningD1)/nD1:0.2%} ({sum(asymmetricTuningD1)}/{nD1})')
print(f'ND1: sym = {sum(symmetricTuningND1)/nND1:0.2%} ({sum(symmetricTuningND1)}/{nND1}), ' + \
      f'asym = {sum(asymmetricTuningND1)/nND1:0.2%} ({sum(asymmetricTuningND1)}/{nND1})')
print(f'FRA types ratio: pValue={pvalFRAtypes:0.4} (oddsRatio={oddsr:0.4})')
print(f'D1neg={nD1neg}  ND1neg={nND1neg}')
print('')

plt.sca(axFRAtypes)
bwidth = 0.7
lineWidth = 1.5
plt.bar(1, 100*nSymmetricTuningD1/nD1, bwidth, ec=figparams.colors['D1'], fc='w', lw=lineWidth)
plt.bar(1, 100*nAsymmetricTuningD1/nD1, bwidth, ec=figparams.colors['D1'], fc=figparams.colors['D1'],
        bottom=100*nSymmetricTuningD1/nD1, lw=lineWidth)
plt.bar(1, 100*nD1neg/nD1, bwidth, ec=figparams.colors['D1'], fc='w',
        bottom=100*(nSymmetricTuningD1+nAsymmetricTuningD1)/nD1, lw=lineWidth)
plt.bar(2, 100*nSymmetricTuningND1/nND1, bwidth, ec=figparams.colors['ND1'], fc='w', lw=lineWidth)
plt.bar(2, 100*nAsymmetricTuningND1/nND1, bwidth, ec=figparams.colors['ND1'], fc=figparams.colors['ND1'],
        bottom=100*nSymmetricTuningND1/nND1, lw=lineWidth)
plt.bar(2, 100*nND1neg/nND1, bwidth, ec=figparams.colors['ND1'], fc='w',
        bottom=100*(nSymmetricTuningND1+nAsymmetricTuningND1)/nND1, lw=lineWidth)
txpos = 2.85
plt.text(txpos, 30, 'Sh', ha='center', va='center', color='k')
plt.text(txpos, 73, 'St', ha='center', va='center', color='k')
plt.text(txpos, 92, 'Su', ha='center', va='center', color='k')
# 2.75: L H N

plt.xticks([1,2])
axFRAtypes.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.25, 2.75])
extraplots.boxoff(axFRAtypes)
plab = plt.ylabel('% low-threshold cells', fontsize=fontSizeLabels, labelpad=0.0)

'''
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

plt.sca(axOSI)
if 1:
    plt.boxplot(osiD1, positions=[1])
    plt.boxplot(osiND1, positions=[2])
else:
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
'''



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
