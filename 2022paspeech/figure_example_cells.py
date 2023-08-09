import os
import sys
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from jaratoolbox import colorpalette as cp
from scipy import signal
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib import patches
from jaratoolbox import extraplots
import studyparams
import figparams
from importlib import reload
reload(figparams)

SAVE_FIGURE = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'figure_example_cells'
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 6]
fontSizeTitles = figparams.fontSizeTitles
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.01, 0.53] # Horiz position for panel labels
labelPosY = [0.975, 0.46]    # Vert position for panel labels

## load database
databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)
FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle=True)

pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
selectivityIndexVOT_FTmax = figData['selectivityIndexVOT_FTmax']
selectivityIndexVOT_FTmin = figData['selectivityIndexVOT_FTmin']
selectivityIndexFT_VOTmax = figData['selectivityIndexFT_VOTmax']
selectivityIndexFT_VOTmin = figData['selectivityIndexFT_VOTmin']
whichFT = figData['whichFT']
whichVOT = figData['whichVOT']
exampleCells = [898, 1155, 162, 687] #VOT selective, FT selective, Mixed, Responsive-nonselective
#exampleCells = [12, 14, 15, 18] #test cells

VOTlabels = ['0', '20', '40', '60']
FTlabels = ['9', '3', '-3', '-9']
colorSounds = cp.TangoPalette['Butter2']
colorsEachVOT = [cp.TangoPalette['Orange1'], cp.TangoPalette['Orange2'],
                 cp.TangoPalette['ScarletRed2'], cp.TangoPalette['ScarletRed3']]
colorsEachFT = [cp.TangoPalette['SkyBlue3'], cp.TangoPalette['SkyBlue1'],
                cp.TangoPalette['Plum1'], cp.TangoPalette['Plum3']]


#timeRange = [-0.3, 0.45]  # In seconds
#timeRange = [-0.15, 0.4]  # In seconds
timeRange = [-0.2, 0.5]  # In seconds
xLims = [-0.15, 0.45]
#xTicks = [0,0.25]
#xTicks = [0, 0.2, 0.4]
xTicks = [0, 0.3]

stimDuration = 0.240 ### THIS SHOULD NOT BE HARDCODED!
colorStim = cp.TangoPalette['Butter3']


def plot_stim(yPos, stimDuration, stimLineWidth=4, stimColor='#edd400'):
    # -- Plot the stimulus --
    #yPos = 1.0*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    pstim = plt.plot([0, stimDuration], 2*[yPos], lw=stimLineWidth, color=stimColor,
                     clip_on=False, solid_capstyle='butt')
    return pstim[0]


plt.clf()

gsMain = gs.GridSpec(2, 2)
gsMain.update(left=0.075, right=0.99, top=0.98, bottom=0.0725, wspace=0.3, hspace=0.3)
#plt.gcf().set_size_inches([14, 12])


for indCell, thisCell in enumerate(exampleCells):
    gsCell = gsMain[indCell].subgridspec(2,2, height_ratios = [0.6, 0.4])
    axRasterVot = plt.subplot(gsCell[0,0])
    axRasterFt = plt.subplot(gsCell[0,1])
    axPsthVot = plt.subplot(gsCell[1,0], sharex = axRasterVot)
    axPsthFt = plt.subplot(gsCell[1,1], sharex = axRasterFt)
    plt.subplots_adjust(wspace = 0.45)

    dbRow = celldb.loc[thisCell]
    oneCell = ephyscore.Cell(dbRow)
    subject = dbRow.subject
    #--Load data FTVOTBorders
    ephysData, bdata = oneCell.load('FTVOTBorders')

    # Align spikes to an event
    spikeTimes = ephysData['spikeTimes']
    eventOnsetTimes = ephysData['events']['stimOn']
    FTParamsEachTrial = bdata['targetFTpercent']


    if (len(FTParamsEachTrial)>len(eventOnsetTimes)) or (len(FTParamsEachTrial)<len(eventOnsetTimes)-1):
        print(f'[{indRow}] Warning! BevahTrials ({len(rateEachTrial)}) and ' +
              f'EphysTrials ({len(eventOnsetTimes)})')
        continue
    if len(FTParamsEachTrial) == len(eventOnsetTimes)-1:
        eventOnsetTimes = eventOnsetTimes[:len(FTParamsEachTrial)]

    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    # Type-sorted rasters -- FTVOTBorders
    VOTParamsEachTrial = bdata['targetVOTpercent']
    possibleVOTParams = np.unique(VOTParamsEachTrial)
    possibleFTParams = np.unique(FTParamsEachTrial)


    trialsEachCond = behavioranalysis.find_trials_each_combination(VOTParamsEachTrial, possibleVOTParams,
                                                                   FTParamsEachTrial, possibleFTParams)
    trialsEachVOT_FTmin = trialsEachCond[:, :, 0]
    trialsEachVOT_FTmax = trialsEachCond[:, :, -1]
    trialsEachFT_VOTmin = trialsEachCond[:, 0, :]
    trialsEachFT_VOTmax = trialsEachCond[:, -1, :]
    nVOT = len(possibleVOTParams)
    nFT = len(possibleFTParams)
    pointSize = 2
    fillWidth = 0.09

    # Raster -- VOT
    plt.sca(axRasterVot)
    if whichVOT[thisCell] == 1: #VOT used for SI == FTmax
        pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                       timeRange, trialsEachVOT_FTmax, colorsEachVOT,
                                                       labels=VOTlabels, fillWidth=fillWidth)
        ymax = np.sum(trialsEachVOT_FTmax)
        #plt.title(rf'$\Delta$VOT, FTmax. SI={np.round(selectivityIndexVOT_FTmax[thisCell], 2)}',
        #          fontsize=fontSizeTitles, fontweight='bold')
    else:
        pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                       timeRange, trialsEachVOT_FTmin, colorsEachVOT,
                                                       labels=VOTlabels, fillWidth=fillWidth)
        #plt.title(rf'$\Delta$VOT, FTmin. SI={np.round(selectivityIndexVOT_FTmin[thisCell], 2)}',
        #           fontsize=fontSizeTitles, fontweight='bold')
        ymax = np.sum(trialsEachVOT_FTmin)
    plt.setp(pRaster, ms=pointSize)
    plt.setp(axRasterVot.get_xticklabels(), visible=False)
    #axRasterVot.set_xticklabels([])  # This would remove them on both shared axes
    plt.ylabel('VOT (ms)', fontsize=fontSizeLabels, fontweight='normal')
    #rect = patches.Rectangle((0, ymax), 0.24, 10, edgecolor='none', facecolor=colorSounds)
    #axRasterVot.add_patch(rect)
    #axRasterVot.text(0, ymax, '', bbox=dict(boxstyle="square", facecolor = colorSounds))

    # Raster -- FT
    plt.sca(axRasterFt)
    if whichFT[thisCell] == 1:
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                      timeRange, trialsEachFT_VOTmax, colorsEachFT,
                                                      labels = FTlabels, fillWidth=fillWidth)
        #plt.title(rf'$\Delta$FT, VOTmax. SI={np.round(selectivityIndexFT_VOTmax[thisCell], 2)}',
        #          fontsize=fontSizeTitles, fontweight='normal')
    else:
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                      timeRange, trialsEachFT_VOTmin, colorsEachFT,
                                                      labels = FTlabels, fillWidth=fillWidth)
        #plt.title(rf'$\Delta$FT, VOTmin. SI={np.round(selectivityIndexFT_VOTmin[thisCell], 2)}',
        #          fontsize=fontSizeTitles, fontweight='normal')
    plt.setp(pRaster, ms=pointSize)
    plt.setp(axRasterFt.get_xticklabels(), visible=False)
    plt.ylabel('FT slope (oct/s)', fontsize=fontSizeLabels, fontweight='normal')


    #-- PSTHs
    lineWidth = 2
    binWidth = 0.010
    timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
    smoothWinSizePsth = 6
    downsampleFactorPsth = 3
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                             timeVec)

    # PSTH -- VOT
    plt.sca(axPsthVot)
    if whichVOT[thisCell] == 1:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachVOT_FTmax,
                             colorsEachVOT, linestyle=None, linewidth=lineWidth)
    else:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachVOT_FTmin,
                             colorsEachVOT, linestyle=None, linewidth=lineWidth)
    plt.ylabel('Firing rate (spk/s)', fontsize=fontSizeLabels, fontweight='normal')
    plt.xlabel('Time (s)', fontsize=fontSizeLabels, fontweight='normal')
    plt.xticks(xTicks)
    extraplots.boxoff(axPsthVot)
    plt.xlim(xLims)
    #PSTHyLims = plt.ylim()
    PSTHyLims = np.ceil(plt.ylim()) * [0,1]
    plt.yticks(PSTHyLims)
    plot_stim(1.1*PSTHyLims[-1], stimDuration)
    plt.ylim(PSTHyLims)


    # PSTH -- FT
    plt.sca(axPsthFt)
    if whichFT[thisCell] == 1:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachFT_VOTmax,
                             colorsEachFT, linestyle=None, linewidth=lineWidth)
    else:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachFT_VOTmin,
                             colorsEachFT, linestyle=None, linewidth=lineWidth)
    plt.xlabel('Time (s)', fontsize=fontSizeLabels, fontweight='normal')
    plt.xticks(xTicks)
    extraplots.boxoff(axPsthFt)
    plt.xlim(xLims)
    #PSTHyLims = np.ceil(plt.ylim()) * [0,1]
    plt.yticks(PSTHyLims)
    plot_stim(1.1*PSTHyLims[-1], stimDuration)
    plt.ylim(PSTHyLims)

    # ***NOTE*** The code above assumes the PSTH for FT is smaller than for VOT (to set ylim).
    #            This may need to change if we use a different set of cells.


plt.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
plt.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
plt.annotate('C', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
plt.annotate('D', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)

#plt.close()
