"""
Figure about responses to amplitude modulated (AM) noise.
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
reload(extraplots)

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_am_responses' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [4.6, 3.8] # In inches
figSize = [7, 4] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

labelPosXtop = [0.01, 0.29, 0.52, 0.75]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.53, 0.75]   # Horiz position for panel labels
labelPosY = [0.94, 0.45]    # Vert position for panel labels

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220213.h5')
celldb = celldatabase.load_hdf(dbPath)

pp = lambda x: f"['{x.subject}', '{x.date}', {x.pdepth}, {x.egroup}, {x.cluster}]"
# pp(celldb.loc[])

selectedD1s = ['d1pi042', '2019-09-11', 3200, 3, 4] # cell 1949
selectedD1t = ['d1pi046', '2020-02-18', 3600, 6, 2] # cell 3094
selectedND1s = ['d1pi047', '2020-02-26', 3500, 4, 3] # cell 3420
selectedND1t = ['d1pi046', '2020-02-18', 3700, 6, 5] # cell 3124 
cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
selectedCells = [selectedD1s, selectedD1t, selectedND1s, selectedND1t]
examplesTypes = ['D1', 'D1', 'ND1', 'ND1']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 1)
gsMain.update(left=0.09, right=0.96, top=0.95, bottom=0.12, hspace=0.6)

gsRast = gsMain[0,0].subgridspec(1, 4, wspace=0.2)
axRastD1s = plt.subplot(gsRast[0, 0])
axRastD1s.annotate('A', xy=(labelPosXtop[0],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axRastD1t = plt.subplot(gsRast[0, 1])
axRastD1t.annotate('B', xy=(labelPosXtop[1],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axRastND1s = plt.subplot(gsRast[0, 2])
axRastND1s.annotate('C', xy=(labelPosXtop[2],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axRastND1t = plt.subplot(gsRast[0, 3])
axRastND1t.annotate('D', xy=(labelPosXtop[3],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axRast = [axRastD1s, axRastD1t, axRastND1s, axRastND1t,]


gsBottom = gsMain[1,0].subgridspec(1, 3, wspace=0.6, width_ratios=[0.6, 0.15, 0.15])
axHist = plt.subplot(gsBottom[0, 0])
axHist.annotate('E', xy=(labelPosXbot[0],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axRS = plt.subplot(gsBottom[0, 1])
axRS.annotate('F', xy=(labelPosXbot[1],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')
axSy = plt.subplot(gsBottom[0, 2])
axSy.annotate('G', xy=(labelPosXbot[2],labelPosY[1]), xycoords='figure fraction',
                fontsize=fontSizePanel, fontweight='bold')


if 1:
    for indType, cellType in enumerate(examplesTypes):
        plt.sca(axRast[indType])
        colorThisType = figparams.colors[cellType]
        from matplotlib import colors as mco
        otherColor = mco.hsv_to_rgb(mco.rgb_to_hsv(mco.to_rgb(colorThisType))*[1,0.4,1]) 
        rasterColors = 6*[colorThisType, otherColor]
        indRow, dbRow = celldatabase.find_cell(celldb, *selectedCells[indType])
        oneCell = ephyscore.Cell(dbRow)

        ephysData, bdata = oneCell.load('am')
        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        if len(bdata['currentFreq']) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(bdata['currentFreq'])]
        timeRange = [-0.2, 0.7]  # In seconds

        rateEachTrial = bdata['currentFreq']
        possibleRate = np.unique(rateEachTrial)
        nRate = len(possibleRate)
        trialsEachCond = behavioranalysis.find_trials_each_type(rateEachTrial, possibleRate)

        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        #yLabels = np.full(nRate, np.nan)
        #yLabels[[0,4,8]] = possibleRate[[0,4,8]]
        yLabels = possibleRate
        pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                       timeRange, trialsEachCond,
                                                       colorEachCond=rasterColors, labels=yLabels)
        plt.setp(pRaster, ms=1)
        plt.xticks([0, 0.5])
        plt.xlabel('Time (s)', fontsize=fontSizeLabels)
        if indType==0:
            plt.ylabel('AM rate (Hz)', fontsize=fontSizeLabels)
            yTicks, lab = plt.yticks()
            rateLabels = [f'{x:0.0f}' for x in possibleRate]
            yTickLabels = nRate*['']
            for lpos in range(0,12,2): #[0,5,10]:
                yTickLabels[lpos] = rateLabels[lpos]
            axRast[indType].set_yticklabels(yTickLabels)
        else:
            axRast[indType].set_yticklabels([])
    plt.show()


# -- Tone response index --
restrictND1 = 0  # True: use only ND1 from tetrodes with D1
toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1 = \
    studyutils.select_cells(celldb, restrictND1=restrictND1)

cellsWithAM = ~np.isnan(celldb['amMinPvalOnset'])
# -- Index for sustained response --
amRespIndex = ( (celldb.amFiringRateBestSustain - celldb.amFiringRateBaseline) /
                (celldb.amFiringRateBestSustain + celldb.amFiringRateBaseline) )

amRespIndexD1 = amRespIndex[indD1 & cellsWithAM]
amRespIndexND1 = amRespIndex[indND1 & cellsWithAM]
amRespIndexEachType = [amRespIndex[indD1], amRespIndex[indND1]]

amRespIndexD1pos = amRespIndexD1[amRespIndexD1>0]
amRespIndexND1pos = amRespIndexND1[amRespIndexND1>0]
amRespIndexD1neg = amRespIndexD1[amRespIndexD1<0]
amRespIndexND1neg = amRespIndexND1[amRespIndexND1<0]
amRespIndexEachTypePos = [amRespIndexD1pos, amRespIndexND1pos]
amRespIndexEachTypeNeg = [amRespIndexD1neg, amRespIndexND1neg]


def line_histogram(data, edges, percent=False, **kwargs):
    hh, ee = np.histogram(data, edges)
    if percent:
        hh = 100 * hh/len(data[~np.isnan(data)])
    hline = plt.step(np.r_[ee, ee[-1]], np.r_[0,hh,0],'-', where='pre', **kwargs)
    return hline
   

plt.sca(axHist)
yLims = [0, 15]
binEdges = np.linspace(-1, 1, 27) #27 32
for indType, cellType in enumerate(cellTypes):
    colorThisType = figparams.colors[cellType]
    hline = line_histogram(amRespIndexEachType[indType], binEdges, lw=2,
                           percent=True, color=colorThisType)
    #hline = line_histogram(np.abs(amRespIndexEachType[indType]), binEdges, lw=2,
    #                       percent=True, color=colorThisType)
    medianPos = np.median(amRespIndexEachTypePos[indType])
    medianNeg = np.median(amRespIndexEachTypeNeg[indType])
    plt.plot(medianPos, 0.95*yLims[-1], 'v', color=colorThisType, mew=2, mfc='none')
    plt.plot(medianNeg, 0.95*yLims[-1], 'v', color=colorThisType, mew=2, mfc='none')
    
plt.ylim(yLims)
#plt.ylim([0,60])
extraplots.boxoff(axHist)
plt.yticks(np.arange(0,20,5))
plt.xlabel('AM sustained response index')
plt.ylabel('% cells')
plt.axvline(0, color='0.9')
#plt.legend(['D1','non-D1'], loc='upper left')
plt.text(-1, 11, 'D1 cells', ha='left', fontweight='bold', fontsize=fontSizeLabels+0,
         color=figparams.colors['D1'])
plt.text(-1, 9, 'non-D1 cells', ha='left', fontweight='bold', fontsize=fontSizeLabels+0,
         color=figparams.colors['ND1'])

plt.show()


uval, pValPos = stats.mannwhitneyu(amRespIndexD1pos, amRespIndexND1pos, alternative='two-sided')
print(f'AM sust index (pos): ' +
      f'D1: {amRespIndexD1pos.median():0.4f} vs ' +
      f'ND1: {amRespIndexND1pos.median():0.4f}   ' +
      f'p = {pValPos:0.6f} (U={uval})')
uval, pValNeg = stats.mannwhitneyu(amRespIndexD1neg, amRespIndexND1neg, alternative='two-sided')
print(f'AM sust index (neg): D1: {amRespIndexD1neg.median():0.4f} vs ' +
      f'ND1: {amRespIndexND1neg.median():0.4f}   ' +
      f'p = {pValNeg:0.6f} (U={uval})')
uval, pValNeg = stats.mannwhitneyu(np.abs(amRespIndexD1[~np.isnan(amRespIndexD1)]),
                                   np.abs(amRespIndexND1[~np.isnan(amRespIndexND1)]),
                                   alternative='two-sided')
print(f'AM sust index (abs): ' +
      f'D1 ({len(amRespIndexD1)}): {np.abs(amRespIndexD1).median():0.4f} vs ' +
      f'ND1 ({len(amRespIndexND1)}): {np.abs(amRespIndexND1).median():0.4f}   ' +
      f'p = {pValNeg:0.6f} (U={uval})')
print()



# -------------------- Rate selectivity ------------------------
maxFiring = celldb['amFiringRateMaxSustain']
minFiring = celldb['amFiringRateMinSustain']
rateSelectivityIndex = (maxFiring - minFiring) / (maxFiring + minFiring)

# NOTE: I need to remove nans because some neurons may be have FR=0 for AM but not tones
amResponsive = amSustainResponsive
rsD1 = rateSelectivityIndex[indD1 & amResponsive & ~np.isnan(rateSelectivityIndex)]
rsND1 = rateSelectivityIndex[indND1 & amResponsive & ~np.isnan(rateSelectivityIndex)]

nD1rs = len(rsD1)
nND1rs = len(rsND1)
medianD1rs = np.median(rsD1)
medianND1rs = np.median(rsND1)

plt.sca(axRS)
markerSize = 3
markerWidth = 0.5
lineWidth = 3
jitter = 0.2
xLine = 0.3*np.array([-1,1])
np.random.seed(0)
xpos = 1 + jitter*(2*np.random.rand(nD1rs)-1)
plt.plot(xpos, rsD1, 'o', mfc='none', mec=figparams.colors['D1'], ms=markerSize, mew=markerWidth)
plt.plot(1+xLine, [medianD1rs, medianD1rs], lw=lineWidth, color=figparams.colors['D1'])
xpos = 2 + jitter*(2*np.random.rand(nND1rs)-1)
plt.plot(xpos, rsND1, 'o', mfc='none', mec=figparams.colors['ND1'], ms=markerSize, mew=markerWidth)
plt.plot(2+xLine, [medianND1rs, medianND1rs], lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axRS.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.5, 2.5])
extraplots.boxoff(axRS)
plt.ylabel('AM rate selectivity', fontsize=fontSizeLabels)

uval, pValRS = stats.mannwhitneyu(rsD1, rsND1, alternative='two-sided')
print(f'AM rate selectivity index: D1 ({nD1rs}): {medianD1rs:0.4f} vs ' +
      f'ND1 ({nND1rs}): {medianND1rs:0.4f}   ' +
      f'p = {pValRS:0.6f} (U={uval})\n')


# -------------------- Max Synchronization Rate ------------------------
maxSyncRate = celldb['amMaxSyncRate']

# NOTE: I need to remove nans because some neurons may be have FR=0 for AM but not tones
amResponsive = amSustainResponsive
syncD1 = maxSyncRate[indD1 & amResponsive & ~np.isnan(maxSyncRate)]
syncND1 = maxSyncRate[indND1 & amResponsive & ~np.isnan(maxSyncRate)]
nD1notSync = sum(syncD1==0)
nND1notSync = sum(syncND1==0)
syncD1 = syncD1[syncD1 > 0]
syncND1 = syncND1[syncND1 > 0]

#logSyncD1 = np.log(syncD1)
#logSyncND1 = np.log(syncND1)

nD1sync = len(syncD1)
nND1sync = len(syncND1)
medianD1sync = np.median(syncD1)
medianND1sync = np.median(syncND1)

yTicks = [4, 8, 16, 32, 64, 128]
#yTicks = np.unique(np.r_[syncD1,syncND1])


plt.sca(axSy)
markerSize = 3
markerWidth = 0.5
lineWidth = 3
jitter = 0.2
spread = 0.07
xLine = 0.3*np.array([-1,1])
np.random.seed(0)
markersD1 = extraplots.spread_plot(1, np.log(syncD1).to_numpy(), spread)
plt.setp(markersD1, ms=markerSize, mew=markerWidth, mfc='none', mec=figparams.colors['D1'])
#plt.plot(1+xLine, np.log([medianD1sync, medianD1sync]), lw=lineWidth, color=figparams.colors['D1'])
markersND1 = extraplots.spread_plot(2, np.log(syncND1).to_numpy(), spread)
plt.setp(markersND1, ms=markerSize, mew=markerWidth, mfc='none', mec=figparams.colors['ND1'])
#plt.plot(2+xLine, np.log([medianND1sync, medianND1sync]), lw=lineWidth, color=figparams.colors['ND1'])
plt.xticks([1,2])
axSy.set_xticklabels(['D1', 'non-D1'], fontsize=fontSizeLabels, rotation=30)
plt.xlim([0.4, 2.5])
plt.yticks(np.log(yTicks))
axSy.set_yticklabels(yTicks)
extraplots.boxoff(axSy)
plt.ylabel('Highest sync. rate (Hz)', fontsize=fontSizeLabels)

uval, pValSync = stats.mannwhitneyu(syncD1, syncND1, alternative='two-sided')
print(f'AM max sync rate: D1 ({nD1sync}): {medianD1sync:0.4f} vs ' +
      f'ND1 ({nND1sync}): {medianND1sync:0.4f}   ' +
      f'p = {pValSync:0.6f} (U={uval})\n')


syncTable = [ [nD1sync, nND1sync], [nD1notSync, nND1notSync] ]
oddsratio, syncPvalFisher = stats.fisher_exact(syncTable, alternative='two-sided')
print(f'Synchronized: \t' +
      f'D1:{nD1sync} ({nD1sync/(nD1sync+nD1notSync):0.2%})  vs  ' +
      f'nD1:{nND1sync} ({nND1sync/(nND1sync+nND1notSync):0.2%})' +
      f'\t p = {syncPvalFisher:0.4f}')


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
