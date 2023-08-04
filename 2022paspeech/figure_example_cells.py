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
figSize = [10, 12]
fontSizeTitles = figparams.fontSizeTitles
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks

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
#exampleCells = [898, 1155, 1054, 687] #VOT selective, FT selective, Mixed, Responsive-nonselective
exampleCells = [12, 14, 15, 18]

VOTlabels = ['0', '20', '40', '60']
FTlabels = ['9', '3', '-3', '-9']
colorSounds = cp.TangoPalette['Butter2']
colorsEachVOT = [cp.TangoPalette['Orange1'], cp.TangoPalette['Orange2'], cp.TangoPalette['ScarletRed2'], cp.TangoPalette['ScarletRed3']]
colorsEachFT = [cp.TangoPalette['SkyBlue3'], cp.TangoPalette['SkyBlue1'], cp.TangoPalette['Plum1'], cp.TangoPalette['Plum3']]
plt.figure()
gsMain = gs.GridSpec(2, 2)
gsMain.update(left=0.075, right=0.925, top=0.925, bottom=0.0725, wspace=0.3, hspace=0.3)
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

    timeRange = [-0.3, 0.45]  # In seconds
    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    # Type-sorted rasters -- FTVOTBorders
    #timeRange = [-0.075, 0.15]
    VOTParamsEachTrial = bdata['targetVOTpercent']
    possibleVOTParams = np.unique(VOTParamsEachTrial)
    possibleFTParams = np.unique(FTParamsEachTrial)


    trialsEachCond = behavioranalysis.find_trials_each_combination(VOTParamsEachTrial, possibleVOTParams, FTParamsEachTrial, possibleFTParams)
    trialsEachVOT_FTmin = trialsEachCond[:, :, 0]
    trialsEachVOT_FTmax = trialsEachCond[:, :, -1]
    trialsEachFT_VOTmin = trialsEachCond[:, 0, :]
    trialsEachFT_VOTmax = trialsEachCond[:, -1, :]
    nVOT = len(possibleVOTParams)
    nFT = len(possibleFTParams)
    pointSize = 4

    # Raster -- VOT
    plt.sca(axRasterVot)
    if whichVOT[thisCell] == 1: #VOT used for SI == FTmax
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmax, colorsEachVOT, labels = VOTlabels)
        ymax = np.sum(trialsEachVOT_FTmax)
        #plt.title(rf'$\Delta$VOT, FTmax. SI={np.round(selectivityIndexVOT_FTmax[thisCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')
    else:
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmin, colorsEachVOT, labels = VOTlabels)
        #plt.title(rf'$\Delta$VOT, FTmin. SI={np.round(selectivityIndexVOT_FTmin[thisCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')
        ymax = np.sum(trialsEachVOT_FTmin)
    plt.setp(pRaster, ms=pointSize)
    plt.xticks([])
    plt.xlim(-0.1,0.35)
    plt.ylabel('VOT (ms)', fontsize=fontSizeLabels, fontweight='bold')
    #rect = patches.Rectangle((0, ymax), 0.24, 10, edgecolor='none', facecolor=colorSounds)
    #axRasterVot.add_patch(rect)
    axRasterVot.text(0, ymax, '', bbox=dict(boxstyle="square", facecolor = colorSounds))



    # Raster -- FT
    plt.sca(axRasterFt)
    if whichFT[thisCell] == 1:
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmax, colorsEachFT, labels = FTlabels)
        #plt.title(rf'$\Delta$FT, VOTmax. SI={np.round(selectivityIndexFT_VOTmax[thisCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')
    else:
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmin, colorsEachFT, labels = FTlabels)
        #plt.title(rf'$\Delta$FT, VOTmin. SI={np.round(selectivityIndexFT_VOTmin[thisCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')
    plt.setp(pRaster, ms=pointSize)
    plt.xticks([])
    plt.xlim(-0.1,0.35)
    plt.ylabel('FT slope (oct/s)', fontsize=fontSizeLabels, fontweight='bold')


    #-- PSTHs
    binWidth = 0.010
    timeRange = [-0.3, 0.45]
    timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
    smoothWinSizePsth = 6
    downsampleFactorPsth = 3
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial, timeVec)
    xTicks = [0,0.25]

    # PSTH -- VOT
    plt.sca(axPsthVot)
    if whichVOT[thisCell] == 1:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachVOT_FTmax, colorsEachVOT, linestyle=None)
    else:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachVOT_FTmin, colorsEachVOT, linestyle=None)
    plt.ylabel('Firing rate (spk/s)', fontsize=fontSizeLabels, fontweight='bold')
    plt.xlabel('Time (s)', fontsize=fontSizeLabels, fontweight='bold')
    plt.xticks(xTicks)

    # PSTH -- FT
    plt.sca(axPsthFt)
    if whichFT[thisCell] == 1:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachFT_VOTmax, colorsEachFT, linestyle=None)
    else:
        extraplots.plot_psth(spikeCountMat/binWidth, smoothWinSizePsth, timeVec, trialsEachFT_VOTmin, colorsEachFT, linestyle=None)
    plt.xlabel('Time (s)', fontsize=fontSizeLabels, fontweight='bold')
    plt.xticks(xTicks)



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)

#plt.close()
