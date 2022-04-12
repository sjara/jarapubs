import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import behavioranalysis
from jaratoolbox import spikesanalysis
from jaratoolbox import celldatabase
from scipy import stats
import pandas as pd
import figparams
reload(figparams)

FIGNAME = 'figure_am'
dataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
# dbPath = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, 'celldatabase_ALLCELLS_MODIFIED_CLU.h5')
dbPath = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, 'celldatabase_calculated_columns.h5')
# dbPath = '/tmp/database_with_pulse_responses.h5'

outputDir='/tmp'

# db = pd.read_hdf(dbPath, key='dataframe')
db = celldatabase.load_hdf(dbPath)
# db = db.query("subject=='pinp015'")
# goodLaser = db.query('pulsePval<0.05 and pulseZscore>0 and trainRatio>0.8')
# goodLaser = db[db['taggedCond']==0]
goodISI = db.query('isiViolations<0.02 or modifiedISI<0.02')
goodShape = goodISI.query('spikeShapeQuality > 2')
goodLaser = goodShape.query("autoTagged==1 and subject != 'pinp018'")
# goodLaser = goodShape.query("autoTagged==1 and subject != 'pinp018' and summaryPulseLatency < 0.01")
goodNSpikes = goodLaser.query('nSpikes>2000')
goodPulseLatency = goodNSpikes.query('summaryPulseLatency<0.006')

# goodSoundResponsiveBool = (~pd.isnull(goodNSpikes['BW10'])) | (~pd.isnull(goodNSpikes['highestSyncCorrected'])) | (goodNSpikes['noiseZscore']<0.05)
# goodSoundResponsive = goodNSpikes[goodSoundResponsiveBool]

dbToUse = goodPulseLatency

ac = dbToUse.groupby('brainArea').get_group('rightAC')
thal = dbToUse.groupby('brainArea').get_group('rightThal')

np.random.seed(1)

messages = []

def jitter(arr, frac):
    jitter = (np.random.random(len(arr))-0.5)*2*frac
    jitteredArr = arr + jitter
    return jitteredArr

def medline(yval, midline, width, color='k', linewidth=3):
    start = midline-(width/2)
    end = midline+(width/2)
    plt.plot([start, end], [yval, yval], color=color, lw=linewidth)

PANELS=[1,1,1,1,1]

SAVE_FIGURE = 1
# outputDir = '/tmp/'
# outputDir = '/mnt/jarahubdata/reports/nick/20171218_all_2018thstr_figures'
# outputDir = figparams.FIGURE_OUTPUT_DIR
figFilename = 'plots_am_tuning' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
# figFormat = 'pdf' # 'pdf' or 'svg'
# figSize = [13,8] # In inches

fullPanelWidthInches = 6.9
figSizeFactor = 2
figWidth = fullPanelWidthInches * (figSizeFactor)
figHeight = figWidth / 1.625
figSize = [figWidth, figHeight] # In inches


thalHistColor = '0.4'
acHistColor = '0.4'

fontSizeLabels = figparams.fontSizeLabels * figSizeFactor
fontSizeTicks = figparams.fontSizeTicks * figSizeFactor
fontSizePanel = figparams.fontSizePanel * figSizeFactor
fontSizeTitles = 12

#Params for extraplots significance stars
fontSizeNS = figparams.fontSizeNS
fontSizeStars = figparams.fontSizeStars
starHeightFactor = figparams.starHeightFactor
starGapFactor = figparams.starGapFactor
starYfactor = figparams.starYfactor
dotEdgeColor = figparams.dotEdgeColor
dataMS = 6

labelPosX = [0.02, 0.35, 0.68, 0.85]   # Horiz position for panel labels
labelPosY = [0.46, 0.96]    # Vert position for panel labels

# Define colors, use figparams
laserColor = figparams.colp['blueLaser']
colorATh = figparams.cp.TangoPalette['SkyBlue2']
colorAC = figparams.cp.TangoPalette['ScarletRed1']

fig = plt.gcf()
plt.clf()
fig.set_facecolor('w')

gs = gridspec.GridSpec(2, 6)
gs.update(left=0.05, right=0.98, top=0.94, bottom=0.10, wspace=0.8, hspace=0.5)

#Load example data
exampleDataPath = os.path.join(dataDir, 'data_am_examples.npz')
exampleData = np.load(exampleDataPath)

exampleFreqEachTrial = exampleData['exampleFreqEachTrial'].item()
exampleSpikeTimes = exampleData['exampleSpikeTimes'].item()
exampleTrialIndexForEachSpike = exampleData['exampleTrialIndexForEachSpike'].item()
exampleIndexLimitsEachTrial = exampleData['exampleIndexLimitsEachTrial'].item()

def plot_example_with_rate(subplotSpec, exampleName, color='k'):
    fig = plt.gcf()

    gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=subplotSpec, wspace=-0.45, hspace=0.0 )

    specRaster = gs[0:2]
    axRaster = plt.Subplot(fig, specRaster)
    fig.add_subplot(axRaster)

    spikeTimes = exampleSpikeTimes[exampleName]
    indexLimitsEachTrial = exampleIndexLimitsEachTrial[exampleName]
    timeRange = [-0.2, 0.7]
    freqEachTrial = exampleFreqEachTrial[exampleName]
    possibleFreq = np.unique(freqEachTrial)
    freqLabels = ['{0:.0f}'.format(freq) for freq in possibleFreq]
    trialsEachCondition = behavioranalysis.find_trials_each_type(freqEachTrial,possibleFreq)
    pRaster, hCond, zline = extraplots.raster_plot(spikeTimes, indexLimitsEachTrial,
                                                   timeRange, trialsEachCondition, labels=freqLabels)
    plt.setp(pRaster, ms=figparams.rasterMS)

    blankLabels = ['']*11
    for labelPos in [0, 5, 10]:
        blankLabels[labelPos] = freqLabels[labelPos]

    axRaster.set_yticklabels(blankLabels)


    ax = plt.gca()
    ax.set_xticks([0, 0.5])
    ax.set_xlabel('Time from\nsound onset (s)', fontsize=fontSizeLabels, labelpad=-1)
    ax.set_ylabel('AM rate (Hz)', fontsize=fontSizeLabels, labelpad=-5)

    # ax.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
    #             fontsize=fontSizePanel, fontweight='bold')

    countRange = [0.1, 0.5]
    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimes,indexLimitsEachTrial,countRange)
    numSpikesInTimeRangeEachTrial = np.squeeze(spikeCountMat)

    numSpikesInTimeRangeEachTrial = np.squeeze(np.diff(indexLimitsEachTrial,
                                                       axis=0))

    if len(numSpikesInTimeRangeEachTrial) == len(freqEachTrial)+1:
        numSpikesInTimeRangeEachTrial = numSpikesInTimeRangeEachTrial[:-1]
    conditionMatShape = np.shape(trialsEachCondition)
    numRepeats = np.product(conditionMatShape[1:])
    nSpikesMat = np.reshape(numSpikesInTimeRangeEachTrial.repeat(numRepeats),
                            conditionMatShape)
    spikesFilteredByTrialType = nSpikesMat * trialsEachCondition
    avgSpikesArray = np.sum(spikesFilteredByTrialType, 0) / np.sum(
        trialsEachCondition, 0).astype('float')/np.diff(np.array(countRange))
    stdSpikesArray = np.std(spikesFilteredByTrialType, 0)/np.diff(np.array(countRange))

    specRate = gs[3]
    axRate = plt.Subplot(fig, specRate)
    fig.add_subplot(axRate)

    nRates = len(possibleFreq)
    plt.hold(True)
    plt.plot(avgSpikesArray, range(nRates), 'ro-', mec='none', ms=6, lw=3, color=color)
    plt.plot(avgSpikesArray-stdSpikesArray, range(len(possibleFreq)), 'k:')
    plt.plot(avgSpikesArray+stdSpikesArray, range(len(possibleFreq)), 'k:')
    axRate.set_ylim([-0.5, nRates-0.5])
    axRate.set_yticks(range(nRates))
    axRate.set_yticklabels([])

    #ax = plt.gca()
    axRate.set_xlabel('Firing rate\n(spk/s)', fontsize = fontSizeLabels, labelpad=-1)
    extraplots.boxoff(axRate)
    # extraplots.boxoff(ax, keep='right')
    return (axRaster, axRate)

spec = gs[0, 0:2]
if PANELS[0]:
    (axRaster, axRate) = plot_example_with_rate(spec, 'Thal0', color=colorATh)
    axRaster.set_title('ATh:Str example 1', fontsize=fontSizeTitles)
    axRate.set_xlim([0,30])
    axRate.set_xticks([0,30])
    extraplots.set_ticks_fontsize(axRate, fontSizeTicks)
    extraplots.set_ticks_fontsize(axRaster, fontSizeTicks)
# ax = plt.gc
axRaster.annotate('A', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

# -- Panel: Thalamus less synchronized --
# axWide = plt.subplot(gs[0, 1])
spec = gs[0, 2:4]
if PANELS[1]:
    (axRaster, axRate) = plot_example_with_rate(spec, 'Thal1', color=colorATh)
    axRaster.set_title('ATh:Str example 2', fontsize=fontSizeTitles)
    axRate.set_xlim([0,35])
    axRate.set_xticks([0,35])
    extraplots.set_ticks_fontsize(axRate, fontSizeTicks)
    extraplots.set_ticks_fontsize(axRaster, fontSizeTicks)
axRaster.annotate('B', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

spec = gs[1, 0:2]
# axSharp.annotate('D', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction',
#              fontsize=fontSizePanel, fontweight='bold')
if PANELS[2]:
    (axRaster, axRate) = plot_example_with_rate(spec, 'AC0', color=colorAC)
    axRaster.set_title('AC:Str example 1', fontsize=fontSizeTitles)
    axRate.set_xlim([0, 12])
    axRate.set_xticks([0, 12])
    extraplots.set_ticks_fontsize(axRate, fontSizeTicks)
    extraplots.set_ticks_fontsize(axRaster, fontSizeTicks)

axRaster.annotate('E', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
# -- Panel: Cortex less synchronized --
# axWide = plt.subplot(gs[1, 1])
spec = gs[1, 2:4]
# axWide.annotate('E', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction',
#              fontsize=fontSizePanel, fontweight='bold')
if PANELS[3]:
    (axRaster, axRate) = plot_example_with_rate(spec, 'AC1', color=colorAC)
    axRaster.set_title('AC:Str example 2', fontsize=fontSizeTitles)
    axRate.set_xlim([0, 12])
    axRate.set_xticks([0, 12])
    extraplots.set_ticks_fontsize(axRate, fontSizeTicks)
    extraplots.set_ticks_fontsize(axRaster, fontSizeTicks)
axRaster.annotate('F', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')



################### Highest Sync #####################


popStatCol = 'highestSyncCorrected'
acPopStat = ac[popStatCol][pd.notnull(ac[popStatCol])]
thalPopStat = thal[popStatCol][pd.notnull(thal[popStatCol])]

acPopStat = acPopStat[acPopStat>0]
thalPopStat = thalPopStat[thalPopStat>0]

# possibleFreqLabels = ["{0:.1f}".format(freq) for freq in np.unique(thalPopStat)]
ytickLabels = [4, 8, 16, 32, 64, 128]
yticks = np.log(ytickLabels)

acPopStat = np.log(acPopStat)
thalPopStat = np.log(thalPopStat)

axSummary = plt.subplot(gs[0, 5])
spacing = 0.07
plt.sca(axSummary)

# pos = jitter(np.ones(len(thalPopStat))*0, 0.20)
# axSummary.plot(pos, thalPopStat, 'o', mec = colorATh, mfc = 'None', alpha=0.5)
plt.hold(1)
markers = extraplots.spread_plot(0, thalPopStat, spacing)
plt.setp(markers, mec = colorATh, mfc = 'None')
plt.setp(markers, ms = dataMS)

plt.hold(1)
medline(np.median(thalPopStat), 0, 0.5)
plt.hold(1)

# pos = jitter(np.ones(len(acPopStat))*1, 0.20)
# axSummary.plot(pos, acPopStat, 'o', mec = colorAC, mfc = 'None', alpha=0.5)
markers = extraplots.spread_plot(1, acPopStat, spacing)
plt.setp(markers, mec = colorAC, mfc = 'None')
plt.setp(markers, ms = dataMS)

plt.hold(1)
medline(np.median(acPopStat), 1, 0.5)
plt.hold(1)

axSummary.set_yticks(yticks)
axSummary.set_yticklabels(ytickLabels)


# tickLabels = ['ATh:Str\nn={}'.format(len(thalPopStat)), 'AC:Str\nn={}'.format(len(acPopStat))]
tickLabels = ['ATh:Str', 'AC:Str']
axSummary.set_xticks(range(2))
axSummary.set_xticklabels(tickLabels, rotation=45)
axSummary.set_xlim([-0.5, 1.5])
extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)
extraplots.boxoff(axSummary)
# axSummary.set_yticks(np.unique(thalPopStat))
# axSummary.set_yticklabels(possibleFreqLabels)
# axSummary.set_ylim([-0.001, 0.161])


yDataMax = max([max(acPopStat), max(thalPopStat)])
yStars = yDataMax + yDataMax*starYfactor
yStarHeight = (yDataMax*starYfactor)*starHeightFactor

zVal, pVal = stats.mannwhitneyu(thalPopStat, acPopStat)
messages.append("{} p={}".format(popStatCol, pVal))
# if pVal < 0.05:
#     extraplots.new_significance_stars([0, 1], np.log(170), np.log(1.1), starMarker='*',
#                                         fontSize=fontSizeStars, gapFactor=starGapFactor)
# else:
#     extraplots.new_significance_stars([0, 1], np.log(170), np.log(1.1), starMarker='n.s.',
#                                         fontSize=fontSizeStars, gapFactor=starGapFactor)
starString = None if pVal<0.05 else 'n.s.'
extraplots.significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
                              starSize=fontSizeStars, starString=starString,
                              gapFactor=starGapFactor)

axSummary.set_ylim([np.log(3.6), np.log(150)])
axSummary.set_ylabel('Highest AM sync. rate (Hz)', labelpad=-1, fontsize=fontSizeLabels)
plt.hold(1)

################### Percent non-sync #####################
# axSummary = plt.subplot(gs[0, 5])

pieChartGS = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=gs[0,4])

axThalPie = plt.subplot(pieChartGS[0, 0])
axACPie = plt.subplot(pieChartGS[1, 0])

annotateX = 0.2
annotateY = np.array([-0.3, -0.45]) + 0.1
rectX = annotateX-0.2
# rectY = [0.5825, 0.55]
rectY = annotateY
rectWidth = 0.15
rectHeight = 0.1

#TODO: Move these to the axis transform
axACPie.annotate('Unsynchronized', xy=(annotateX, annotateY[0]), xycoords='axes fraction', fontsize=fontSizeTicks-2)
axACPie.annotate('Synchronized', xy=(annotateX, annotateY[1]), xycoords='axes fraction', fontsize=fontSizeTicks-2)

fig = plt.gcf()
rect1 = mpatches.Rectangle(xy=(rectX, rectY[0]), width=rectWidth, height=rectHeight, fc = 'w', ec='k', clip_on=False,
                           transform=axACPie.transAxes)
rect2 = mpatches.Rectangle(xy=(rectX, rectY[1]), width=rectWidth, height=rectHeight, fc = 'k', ec='k', clip_on=False,
                           transform=axACPie.transAxes)

axACPie.add_patch(rect1)
axACPie.add_patch(rect2)


popStatCol = 'highestSyncCorrected'
acPopStat = ac[popStatCol][pd.notnull(ac[popStatCol])]
acPopStat = acPopStat[pd.notnull(acPopStat)]
thalPopStat = thal[popStatCol][pd.notnull(thal[popStatCol])]
thalPopStat = thalPopStat[pd.notnull(thalPopStat)]

acSyncN = len(acPopStat[acPopStat > 0])
acNonSyncN = len(acPopStat[acPopStat == 0])
acSyncFrac = acSyncN/float(acSyncN + acNonSyncN)
acNonSyncFrac = acNonSyncN/float(acSyncN + acNonSyncN)

pieWedges = axACPie.pie([acNonSyncFrac, acSyncFrac], colors=['w', colorAC], shadow=False, startangle=0)
for wedge in pieWedges[0]:
    wedge.set_edgecolor(colorAC)

# axACPie.annotate('Non-Sync\n{}%'.format(int(100*acNonSyncFrac)), xy=[0.8, 0.8], rotation=0, fontweight='bold', textcoords='axes fraction')
# axACPie.annotate('Sync\n{}%'.format(int(100*acSyncFrac)), xy=[-0.05, -0.05], rotation=0, fontweight='bold', textcoords='axes fraction')
fontSizePercent = 12
axACPie.annotate('{:0.0f}%'.format(np.round(100*acNonSyncFrac)), xy=[0.48, 0.6], rotation=0,
                 fontweight='regular', textcoords='axes fraction', fontsize=fontSizePercent)
axACPie.annotate('{:0.0f}%'.format(np.round(100*acSyncFrac)), xy=[0.25, 0.25], rotation=0,
                 fontweight='bold', textcoords='axes fraction', fontsize=fontSizePercent, color='w')
axACPie.set_aspect('equal')

thalSyncN = len(thalPopStat[thalPopStat > 0])
thalNonSyncN = len(thalPopStat[thalPopStat == 0])
thalSyncFrac = thalSyncN/float(thalSyncN + thalNonSyncN)
thalNonSyncFrac = thalNonSyncN/float(thalSyncN + thalNonSyncN)

pieWedges = axThalPie.pie([thalNonSyncFrac, thalSyncFrac], colors=['w', colorATh], shadow=False, startangle=0)
for wedge in pieWedges[0]:
    wedge.set_edgecolor(colorATh)

# axThalPie.annotate('Non-Sync\n{}%'.format(int(100*thalNonSyncFrac)), xy=[0.8, 0.8], rotation=0, fontweight='bold', textcoords='axes fraction')
# axThalPie.annotate('Sync\n{}%'.format(int(100*thalSyncFrac)), xy=[-0.05, -0.05], rotation=0, fontweight='bold', textcoords='axes fraction')
axThalPie.annotate('{:0.0f}%'.format(np.round(100*thalNonSyncFrac)), xy=[0.57, 0.525], rotation=0,
                 fontweight='regular', textcoords='axes fraction', fontsize=fontSizePercent)
axThalPie.annotate('{:0.0f}%'.format(np.round(100*thalSyncFrac)), xy=[0.2, 0.3], rotation=0,
                   fontweight='bold', textcoords='axes fraction', fontsize=fontSizePercent, color='w')
axThalPie.set_aspect('equal')

oddsratio, pValue = stats.fisher_exact([[acSyncN, thalSyncN],
                                        [acNonSyncN, thalNonSyncN]])
print "AC: {} Nonsync / {} total".format(acNonSyncN, acSyncN+acNonSyncN)
print "Thal: {} Nonsync / {} total".format(thalNonSyncN, thalSyncN+thalNonSyncN)
print "p-Val for fisher exact test: {}".format(pValue)
if pValue < 0.05:
    starMarker = '*'
else:
    starMarker = 'n.s.'


axThalPie.annotate('C', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
axThalPie.annotate('D', xy=(labelPosX[3],labelPosY[1]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
axThalPie.annotate('G', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')
axThalPie.annotate('H', xy=(labelPosX[3],labelPosY[0]), xycoords='figure fraction',
             fontsize=fontSizePanel, fontweight='bold')

xBar = -2
#FarUntagged, CloseUntagged, tagged
yCircleCenters = [0, 3]
xTickWidth = 0.2
yGapWidth = 0.5

##### -- Plotting stars for the pie charts -- ######
# def plot_y_lines_with_ticks(ax, x, y1, y2, gapwidth, tickwidth, color='k', starMarker="*", fontSize=9):
#     ax.plot([x, x], [y1, np.mean([y1, y2])-(gapwidth/2)], '-', clip_on=False, color=color)
#     ax.hold(1)
#     ax.plot([x, x], [np.mean([y1, y2])+(gapwidth/2), y2], '-', clip_on=False, color=color)
#     ax.plot([x, x+xTickWidth], [y1, y1], '-', clip_on=False, color=color)
#     ax.plot([x, x+xTickWidth], [y2, y2], '-', clip_on=False, color=color)

#     if starMarker == '*':
#         ax.plot(x, np.mean([y1, y2]), starMarker, clip_on=False, ms=fontSizeStars, mfc='k', mec='None')
#     else:
#         ax.text(x, np.mean([y1, y2]), starMarker, fontsize=fontSize, rotation='vertical', va='center', ha='center', clip_on=False)


# plot_y_lines_with_ticks(axACPie, xBar, yCircleCenters[0], yCircleCenters[1],
#                         yGapWidth, xTickWidth, starMarker=starMarker)

#####################################################

# width = 0.5
# plt.hold(1)
# loc = [1, 2]
# # axSummary.bar(loc[0]-width/2, thalNonSyncPercent, width, color=colorATh)
# # axSummary.bar(loc[0]-width/2, thalSyncPercent, width, bottom=thalNonSyncPercent, color=colorATh, alpha=0.5)
# # axSummary.bar(loc[1]-width/2, acNonSyncPercent, width, color=colorAC)
# # axSummary.bar(loc[1]-width/2, acSyncPercent, width, bottom=acNonSyncPercent, color=colorAC, alpha=0.5)
# axSummary.bar(loc[0], thalNonSyncPercent, width, color=colorATh)
# axSummary.bar(loc[0], thalSyncPercent, width, bottom=thalNonSyncPercent, color=colorATh, alpha=0.5)
# axSummary.bar(loc[1], acNonSyncPercent, width, color=colorAC)
# axSummary.bar(loc[1], acSyncPercent, width, bottom=acNonSyncPercent, color=colorAC, alpha=0.5)
# extraplots.boxoff(axSummary)

# extraplots.new_significance_stars([1, 2], 105, 2.5, starMarker='*',
#                                     fontSize=fontSizeStars, gapFactor=starGapFactor)

# axSummary.text(2.65, 30, 'Non-Sync.', rotation=90, fontweight='bold')
# axSummary.text(2.65, 75, 'Sync.', rotation=90, fontweight='bold', color='0.5')

# axSummary.set_xlim([0.5, 2.6])
# # extraplots.boxoff(axSummary)
# axSummary.set_ylim([0, 100.5])
# axSummary.set_xticks([1, 2])
# tickLabels = ['ATh:Str', 'AC:Str']
# axSummary.set_xticklabels(tickLabels)
# axSummary.set_ylabel('% neurons', labelpad=-5)

##########################################################


########## Discrimination of Rate #############
# dbPathRate = os.path.join(dataDir, 'celldatabase_with_am_discrimination_accuracy.h5')
dbPathRate = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, 'celldatabase_calculated_columns.h5')
# dataframeRate = pd.read_hdf(dbPathRate, key='dataframe')
dataframeRate = celldatabase.load_hdf(dbPathRate)
#Show population distributions
colorATh = figparams.cp.TangoPalette['SkyBlue2']
colorAC = figparams.cp.TangoPalette['ScarletRed1']
# dataMS=5

goodISIRate = dataframeRate.query('isiViolations<0.02 or modifiedISI<0.02')
goodShapeRate = goodISIRate.query('spikeShapeQuality > 2')
goodLaserRate = goodShapeRate.query("autoTagged==1 and subject != 'pinp018'")
goodNSpikesRate = goodLaserRate.query('nSpikes>2000')

goodPulseLatency = goodNSpikesRate.query('summaryPulseLatency<0.01')

dbToUse = goodPulseLatency

acRate = dbToUse.groupby('brainArea').get_group('rightAC')
thalRate = dbToUse.groupby('brainArea').get_group('rightThal')

popStatCol = 'rateDiscrimAccuracy'
# popStatCol = 'accuracySustained'
acPopStat = acRate[popStatCol][pd.notnull(acRate[popStatCol])]
thalPopStat = thalRate[popStatCol][pd.notnull(thalRate[popStatCol])]

# plt.clf()
# axSummary = plt.subplot(111)
axSummary = plt.subplot(gs[1, 4])

jitterFrac = 0.2
pos = jitter(np.ones(len(thalPopStat))*0, jitterFrac)
axSummary.plot(pos, thalPopStat, 'o', mec = colorATh, mfc = 'None', alpha=1, ms=dataMS)
medline(np.median(thalPopStat), 0, 0.5)
pos = jitter(np.ones(len(acPopStat))*1, jitterFrac)
axSummary.plot(pos, acPopStat, 'o', mec = colorAC, mfc = 'None', alpha=1, ms=dataMS)
medline(np.median(acPopStat), 1, 0.5)
tickLabels = ['ATh:Str'.format(len(thalPopStat)), 'AC:Str'.format(len(acPopStat))]
axSummary.set_xticks(range(2))
axSummary.set_xticklabels(tickLabels, rotation=45)
extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)
axSummary.set_ylim([0.5, 1])
yticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1]
axSummary.set_yticks(yticks)
ytickLabels = ['50', '', '', '', '', '100']
axSummary.set_yticklabels(ytickLabels)
axSummary.set_ylabel('Discrimination accuracy\nof AM rate (%)', fontsize=fontSizeLabels, labelpad=-12)
# extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)

zstat, pVal = stats.mannwhitneyu(thalPopStat, acPopStat)

messages.append("{} p={}".format("Rate discrimination accuracy", pVal))
messages.append("{} ATh n={}, AC n={}".format(popStatCol, len(thalPopStat), len(acPopStat)))

# plt.title('p = {}'.format(np.round(pVal, decimals=5)))

# axSummary.annotate('C', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
#             fontsize=fontSizePanel, fontweight='bold')

# starHeightFactor = 0.2
# starGapFactor = 0.3
# starYfactor = 0.1
yDataMax = max([max(acPopStat), max(thalPopStat)]) - 0.025
yStars = yDataMax + yDataMax*starYfactor
yStarHeight = (yDataMax*starYfactor)*starHeightFactor

starString = None if pVal<0.05 else 'n.s.'
# fontSizeStars = 9
extraplots.significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
                            starSize=fontSizeStars, starString=starString,
                            gapFactor=starGapFactor)
extraplots.boxoff(axSummary)
plt.hold(1)


########## Discrimination of Phase #############
dbPathPhase = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, 'celldatabase_calculated_columns.h5')
# dbPhase = pd.read_hdf(dbPathPhase, key='dataframe')
dbPhase = celldatabase.load_hdf(dbPathPhase)

goodISIPhase = dbPhase.query('isiViolations<0.02 or modifiedISI<0.02')
goodShapePhase = goodISIPhase.query('spikeShapeQuality > 2')
goodLaserPhase = goodShapePhase.query("autoTagged==1 and subject != 'pinp018'")
goodNSpikesPhase = goodLaserPhase.query('nSpikes>2000')

acPhase = goodNSpikesPhase.groupby('brainArea').get_group('rightAC')
thalPhase = goodNSpikesPhase.groupby('brainArea').get_group('rightThal')

possibleRateKeys = np.array([4, 5, 8, 11, 16, 22, 32, 45, 64, 90, 128])
ratesToUse = possibleRateKeys
keys = ['phaseDiscrimAccuracy_{}Hz'.format(rate) for rate in ratesToUse]

acData = np.full((len(acPhase), len(ratesToUse)), np.nan)
thalData = np.full((len(thalPhase), len(ratesToUse)), np.nan)

for externalInd, (indRow, row) in enumerate(acPhase.iterrows()):
    for indKey, key in enumerate(keys):
        acData[externalInd, indKey] = row[key]

for externalInd, (indRow, row) in enumerate(thalPhase.iterrows()):
    for indKey, key in enumerate(keys):
        thalData[externalInd, indKey] = row[key]

acMeanPerCell = np.nanmean(acData, axis=1)
acMeanPerCell = acMeanPerCell[~np.isnan(acMeanPerCell)]
thalMeanPerCell = np.nanmean(thalData, axis=1)
thalMeanPerCell = thalMeanPerCell[~np.isnan(thalMeanPerCell)]

# plt.clf()

axSummary = plt.subplot(gs[1, 5])

jitterFrac = 0.2
pos = jitter(np.ones(len(thalMeanPerCell))*0, jitterFrac)
axSummary.plot(pos, thalMeanPerCell, 'o', mec = colorATh, mfc = 'None', alpha=1, ms=dataMS)
medline(np.median(thalMeanPerCell), 0, 0.5)
pos = jitter(np.ones(len(acMeanPerCell))*1, jitterFrac)
axSummary.plot(pos, acMeanPerCell, 'o', mec = colorAC, mfc = 'None', alpha=1, ms=dataMS)
medline(np.median(acMeanPerCell), 1, 0.5)
tickLabels = ['ATh:Str'.format(len(thalMeanPerCell)), 'AC:Str'.format(len(acMeanPerCell))]
axSummary.set_xticks(range(2))
axSummary.set_xticklabels(tickLabels, rotation=45)
extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)
axSummary.set_ylim([0.5, 1])
yticks = [0.5, 0.6, 0.7, 0.8, 0.9, 1]
axSummary.set_yticks(yticks)
ytickLabels = ['50', '', '', '', '', '100']
axSummary.set_yticklabels(ytickLabels)
# axSummary.set_yticklabels(map(str, [50, 60, 70, 80, 90, 100]))
axSummary.set_ylabel('Discrimination accuracy\nof AM phase (%)', fontsize=fontSizeLabels, labelpad=-12)


zstat, pVal = stats.mannwhitneyu(thalMeanPerCell, acMeanPerCell)

messages.append("{} p={}".format("Phase discrimination accuracy", pVal))
messages.append("{} ATh n={}, AC n={}".format("Phase discrimination accuracy", len(thalPopStat), len(acPopStat)))

# plt.title('p = {}'.format(np.round(pVal, decimals=5)))

# axSummary.annotate('C', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
#             fontsize=fontSizePanel, fontweight='bold')

# starHeightFactor = 0.2
# starGapFactor = 0.3
# starYfactor = 0.1
yDataMax = max([max(acMeanPerCell), max(thalMeanPerCell)])
yStars = yDataMax + yDataMax*starYfactor
yStarHeight = (yDataMax*starYfactor)*starHeightFactor

starString = None if pVal<0.05 else 'n.s.'
# fontSizeStars = 9
extraplots.significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
                            starSize=fontSizeStars, starString=starString,
                            gapFactor=starGapFactor)


extraplots.boxoff(axSummary)







plt.show()
print "\nSTATISTICS:\n"
for message in messages:
    print(message)
print "\n"

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)

'''
################### Mutual info PHASE #####################
# popStatCol = 'mutualInfoPerSpike'
# acPopStat = ac[popStatCol][pd.notnull(ac[popStatCol])]
# thalPopStat = thal[popStatCol][pd.notnull(thal[popStatCol])]

axSummary = plt.subplot(gs[1, 5])

possibleFreqKeys = [4, 5, 8, 11, 16, 22, 32, 45, 64, 90, 128]

# dataframe = dataframe.query("pulsePval<0.05 and trainRatio>0.8")
# ac = dataframe.groupby('brainArea').get_group('rightAC')
# thal = dataframe.groupby('brainArea').get_group('rightThal')

keys = ['mutualInfoPhase_{}Hz'.format(rate) for rate in possibleFreqKeys]

acData = np.full((len(ac), len(possibleFreqKeys)), np.nan)
thalData = np.full((len(thal), len(possibleFreqKeys)), np.nan)

for externalInd, (indRow, row) in enumerate(ac.iterrows()):
    for indKey, key in enumerate(keys):
        acData[externalInd, indKey] = row[key]

for externalInd, (indRow, row) in enumerate(thal.iterrows()):
    for indKey, key in enumerate(keys):
        thalData[externalInd, indKey] = row[key]

acData[acData<0]=0
thalData[thalData<0]=0

allPval = []
for indCol, freqKey in enumerate(possibleFreqKeys):
    acDataThisFreq = acData[:,indCol][np.logical_not(np.isnan(acData[:,indCol]))]
    thalDataThisFreq = thalData[:,indCol][np.logical_not(np.isnan(thalData[:,indCol]))]
    zStat, pVal = stats.ranksums(acDataThisFreq, thalDataThisFreq)
    allPval.append(int(pVal<0.05))
    print "{}Hz, p={}".format(freqKey, pVal)

acMean = np.nanmean(acData, axis=0)
# acMean = np.nanmedian(acData, axis=0)
acStd = np.nanstd(acData, axis=0)

thalMean = np.nanmean(thalData, axis=0)
# thalMean = np.nanmedian(thalData, axis=0)
thalStd = np.nanstd(thalData, axis = 0)

axSummary.plot(acMean, '-', color=colorAC, label='AC:Str')
# plt.fill_between(range(len(possibleFreqKeys)), acMean+acStd/numAC, acMean-acStd/numAC, color='r', alpha=0.5)
plt.hold(1)
axSummary.plot(thalMean, '-', color=colorATh, label="ATh:Str")
# plt.fill_between(range(len(possibleFreqKeys)), thalMean+thalStd/numThal, thalMean-thalStd/numThal, color='b', alpha=0.5)
axSummary.set_xticks(range(len(possibleFreqKeys))[::2])
axSummary.set_xticklabels(possibleFreqKeys[::2])
axSummary.set_xlabel('AM rate (Hz)')

for indRate, significant in enumerate(allPval):
    if significant:
        axSummary.plot(indRate, np.mean([thalMean[indRate],acMean[indRate]]), "k*")

axSummary.set_ylabel('MI, Spike rate vs. stimulus phase (bits)')
extraplots.boxoff(axSummary)
axSummary.legend()




'''

# ############## Mutual info AM RATE #######################
# popStatCol = 'mutualInfoBCBits'
# # popStatCol = 'mutualInfoPerSpikeBits'
# acPopStat = ac[popStatCol][pd.notnull(ac[popStatCol])]
# thalPopStat = thal[popStatCol][pd.notnull(thal[popStatCol])]

# axSummary = plt.subplot(gs[1, 4])

# #REMOVE NEGATIVE MI VALUES, REPLACE WITH 0
# acPopStat[acPopStat < 0] = 0
# thalPopStat[thalPopStat < 0] = 0
# pos = jitter(np.ones(len(thalPopStat))*0, 0.20)
# axSummary.plot(pos, thalPopStat, 'o', mec = colorATh, mfc = 'None', alpha=1, ms=dataMS)
# medline(np.median(thalPopStat), 0, 0.5)
# pos = jitter(np.ones(len(acPopStat))*1, 0.20)
# axSummary.plot(pos, acPopStat, 'o', mec = colorAC, mfc = 'None', alpha=1, ms=dataMS)
# medline(np.median(acPopStat), 1, 0.5)
# if popStatCol == 'mutualInfoPerSpikeBits':
#     plt.ylabel('MI(Firing rate; AM rate) (bits/spike)', fontsize=fontSizeLabels)
# elif popStatCol == 'mutualInfoBCBits':
#     plt.ylabel('MI(Firing rate; AM rate)\n(bits)', fontsize=fontSizeLabels, labelpad=-10)
# # tickLabels = ['ATh:Str', 'AC:Str']
# # tickLabels = ['ATh:Str\nn={}'.format(len(thalPopStat)), 'AC:Str\nn={}'.format(len(acPopStat))]
# tickLabels = ['ATh:Str'.format(len(thalPopStat)), 'AC:Str'.format(len(acPopStat))]
# axSummary.set_xticks(range(2))
# axSummary.set_xticklabels(tickLabels, rotation=45)
# extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)
# axSummary.set_xlim([-0.5, 1.5])
# extraplots.boxoff(axSummary)

# axSummary.set_yticks([0, 0.25])
# axSummary.set_yticklabels(['0', '0.25'], rotation=90, va='center')
# extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)

# axSummary.set_ylim([-0.001, 0.25])

# zstat, pVal = stats.ranksums(thalPopStat, acPopStat)

# messages.append("{} p={}".format(popStatCol, pVal))

# # plt.title('p = {}'.format(np.round(pVal, decimals=5)))

# axSummary.annotate('C', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
#             fontsize=fontSizePanel, fontweight='bold')

# yDataMax = max([max(acPopStat), max(thalPopStat)])
# yStars = yDataMax + yDataMax*starYfactor
# yStarHeight = (yDataMax*starYfactor)*starHeightFactor

# starString = None if pVal<0.05 else 'n.s.'
# extraplots.significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
#                               starSize=fontSizeStars, starString=starString,
#                               gapFactor=starGapFactor)
# plt.hold(1)


############################################################


# ################### Mutual info PHASE #####################
# # popStatCol = 'mutualInfoPerSpike'
# # acPopStat = ac[popStatCol][pd.notnull(ac[popStatCol])]
# # thalPopStat = thal[popStatCol][pd.notnull(thal[popStatCol])]

# axSummary = plt.subplot(gs[1, 5])

# possibleRateKeys = np.array([4, 5, 8, 11, 16, 22, 32, 45, 64, 90, 128])
# rateThreshold = 22
# ratesToUse = possibleRateKeys[possibleRateKeys>rateThreshold]


# # dataframe = dataframe.query("pulsePval<0.05 and trainRatio>0.8")
# # ac = dataframe.groupby('brainArea').get_group('rightAC')
# # thal = dataframe.groupby('brainArea').get_group('rightThal')

# keys = ['mutualInfoPhase_{}Hz'.format(rate) for rate in ratesToUse]

# acData = np.full((len(ac), len(ratesToUse)), np.nan)
# thalData = np.full((len(thal), len(ratesToUse)), np.nan)

# for externalInd, (indRow, row) in enumerate(ac.iterrows()):
#     for indKey, key in enumerate(keys):
#         acData[externalInd, indKey] = row[key]

# for externalInd, (indRow, row) in enumerate(thal.iterrows()):
#     for indKey, key in enumerate(keys):
#         thalData[externalInd, indKey] = row[key]

# acData = np.nanmean(acData, axis=1)
# thalData = np.nanmean(thalData, axis=1)

# acData[acData<0]=0
# thalData[thalData<0]=0

# thalPopStat = thalData[~np.isnan(thalData)]
# pos = jitter(np.ones(len(thalPopStat))*0, 0.20)
# axSummary.plot(pos, thalPopStat, 'o', mec = colorATh, mfc = 'None', alpha=1, ms=dataMS)
# medline(np.median(thalPopStat), 0, 0.5)

# acPopStat = acData[~np.isnan(acData)]
# pos = jitter(np.ones(len(acPopStat))*1, 0.20)
# axSummary.plot(pos, acPopStat, 'o', mec = colorAC, mfc = 'None', alpha=1, ms=dataMS)
# medline(np.median(acPopStat), 1, 0.5)

# # tickLabels = ['ATh:Str', 'AC:Str']
# # tickLabels = ['ATh:Str\nn={}'.format(len(thalPopStat)), 'AC:Str\nn={}'.format(len(acPopStat))]
# tickLabels = ['ATh:Str', 'AC:Str']
# axSummary.set_xticks(range(2))
# axSummary.set_xticklabels(tickLabels, rotation=45)
# extraplots.set_ticks_fontsize(axSummary, fontSizeLabels)
# axSummary.set_xlim([-0.5, 1.5])
# extraplots.boxoff(axSummary)

# axSummary.set_yticks([0, 0.010])
# axSummary.set_yticklabels(['0', '0.01'], rotation=90, va='center')


# # axSummary.set_ylim([-0.0001, 0.025])

# zstat, pVal = stats.ranksums(thalPopStat, acPopStat)

# messages.append("mutualInfoPhase for rates>22Hz,  p={}".format(pVal))

# # plt.title('p = {}'.format(np.round(pVal, decimals=5)))

# # axSummary.annotate('E', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
# #             fontsize=fontSizePanel, fontweight='bold')

# yDataMax = max([max(acPopStat), max(thalPopStat)])
# # yDataMax = 0.023
# yStars = yDataMax + yDataMax*starYfactor
# yStarHeight = (yDataMax*starYfactor)*starHeightFactor

# starString = None if pVal<0.05 else 'n.s.'
# extraplots.significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
#                               starSize=fontSizeStars, starString=starString,
#                               gapFactor=starGapFactor)
# '''
# if pVal < 0.05:
#     extraplots.new_significance_stars([0, 1], yStars, yStarHeight, starMarker='*',
#                                         fontSize=fontSizeStars, gapFactor=starGapFactor)
# else:
#     extraplots.new_significance_stars([0, 1], yStars, yStarHeight, starMarker='n.s.',
#                                         fontSize=fontSizeStars, gapFactor=starGapFactor)
# '''

# axSummary.set_ylabel('MI(Firing rate; phase)\n(bits)', fontsize=fontSizeLabels, labelpad=-10)
# extraplots.boxoff(axSummary)
