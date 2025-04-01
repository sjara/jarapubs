"""
Compare oddball enhancement under different reagents.
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)
importlib.reload(extraplots)


SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_oddball_enhancement_bars' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [14, 5.5] # In inches

if len(sys.argv)>1:
    trialSubset = sys.argv[1]
    figFilename = figFilename + '_' + trialSubset
    figFormat = 'pdf'
    if trialSubset=='all':
        trialSubset = ''
else:
    trialSubset = ''
if trialSubset not in ['', 'all', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

if len(sys.argv)>2:
    evokedPeriod = sys.argv[2]
    figFilename = figFilename + '_' + evokedPeriod
else:
    evokedPeriod = ''
if evokedPeriod not in ['', 'onset', 'sustained', 'offset']:
    raise ValueError("evokedPeriod must be '', 'onset', 'sustained', 'offset'")


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.01, 0.215, 0.467, 0.705]   # Horiz position for panel labels
labelPosY = [0.95, 0.45]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorsRasterDark = figparams.colors
colorsRasterLight = figparams.colorsLight
pColor = '0.5'
rasterMarkerSize = 0.5

# -- Information about examples --
'''
cellDictChord = {'subject' : 'acid006', 'date' : '2023-03-22',
              'pdepth' : 3000, 'egroup' : 0, 'cluster' : 284} # Cell 86
cellDictFM = {'subject' : 'acid006', 'date' : '2023-03-22',
              'pdepth' : 3000, 'egroup' : 0, 'cluster' : 135} # Cell 46
cellsToPlot = [cellDictFM,  cellDictFM]
'''
'''
cellsToPlot = [133,  142]  #  High:84, 86  Low: 1065
oddballToPlot = ['HighFreq', 'HighFreq']
cellsToPlot = [46,  86]  #  High:84, 86  Low: 1065
oddballToPlot = ['FM_Up', 'HighFreq']
yMaxEachCell = [180, 80]
'''

CASE = 0
if CASE==0:
    cellsToPlot = [86,  46]  #  High:84, 86  Low: 1065
    oddballToPlot = ['HighFreq', 'FM_Up']
    yMaxEachCell = [80, 180]
elif CASE==1:
    cellsToPlot = [403,  46]  #  High:84, 86  Low: 1065  (75:High,)
    oddballToPlot = ['HighFreq', 'FM_Up']
    yMaxEachCell = [110, 180]
elif CASE==2:
    cellsToPlot = [511,  549] 
    oddballToPlot = ['HighFreq', 'FM_Up']
    yMaxEachCell = [50, 50]

    
# -- Raster and PSTH parameters --
timeRange = [-0.3, 0.45]
#timeRange = [-0.32, 0.97]
#timeRangeToShow = [-0.3, 0.45]
binWidth = 0.010
timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
smoothWinSizePsth = 2 
lwPsth = 2.5
downsampleFactorPsth = 1

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldbAll = celldatabase.load_hdf(dbFilename)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball_running.h5')
    celldb = celldatabase.load_hdf(dbFilename)
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball_notrunning.h5')
    celldb = celldatabase.load_hdf(dbFilename)
else:
    celldb = celldbAll
if evokedPeriod != '':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball_{evokedPeriod}.h5')
    celldb = celldatabase.load_hdf(dbFilename)

# -- Calculate oddball enhancement index (OEI) --
#responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=studyparams.FR_THRESHOLD)
for inds, stim in enumerate(studyparams.STIMULI):
    cstim = stim.capitalize()
    oeindex = np.empty((len(celldb), len(studyparams.REAGENTS)))
    for indr, reagent in enumerate(studyparams.REAGENTS):
        standardEvokedFR = celldb[reagent+cstim+'StandardEvokedFiringRate']
        oddballEvokedFR = celldb[reagent+cstim+'OddballEvokedFiringRate']
        celldb[reagent+cstim+'OEI'] = studyutils.modulation_index(oddballEvokedFR, standardEvokedFR)
        #oeindex[:,indr] = studyutils.modulation_index(oddballEvokedFR, standardEvokedFR)

    
def plot_stim(yLims, stimDuration, stimLineWidth=6, stimColor=figparams.colorStim):
    yPos = 0.9*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    pstim = plt.plot([0, stimDuration], 2*[yPos], lw=stimLineWidth, color=stimColor,
                     clip_on=False, solid_capstyle='butt')
    return pstim[0]

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

stimClasses = ['FM', 'pureTone']

# -- Main gridspec --
#gsMain = gridspec.GridSpec(len(stimClasses), 3, width_ratios=[0.16, 0.52, 0.32])
gsMain = gridspec.GridSpec(len(stimClasses), 3)
gsMain.update(left=0.07, right=0.99, top=0.94, bottom=0.08, wspace=0.25, hspace=0.4)

# -- Show panel labels --
# for indp, plabel in enumerate(['A','B','C','D']):
#     plt.figtext(labelPosX[indp], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')
# for indp, plabel in enumerate(['E','F','G','H']):
#     plt.figtext(labelPosX[indp], labelPosY[1], plabel, fontsize=fontSizePanel, fontweight='bold')

panelEachStimType = [[0,3], [0,4], [1,3], [1,4]]

#ax1.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
#             fontsize=fontSizePanel, fontweight='bold')
#axs = np.empty((len(stimClasses), nCols))

if 0:
    # -- Plot examples --
    for indcell, cellInd in enumerate(cellsToPlot):
        gsExamples = gsMain[indcell, 1].subgridspec(1, 2, wspace=0.15)
        #cellInd, dbRow = celldatabase.find_cell(celldb, **cellsToPlot[indcell])
        dbRow = celldb.iloc[cellsToPlot[indcell]]
        oneCell = ephyscore.Cell(dbRow)
        reagentsToPlot = ['saline', 'doi']
        reagentLabels = ['Saline', 'DOI']
        oddballToLoad = oddballToPlot[indcell]
        for indr, reagent in enumerate(reagentsToPlot):
            colorEachCond = figparams.colorsLightDark[reagent]
            (stfeo, ilet, tec, clabels, bdataList) = \
                studyutils.load_oddball_data_one_cell(oneCell, oddballToLoad,
                                                      reagent, timeRange, nPreOdd=2)
            rasterLabels = [s.capitalize() for s in clabels]
            #gsRasterPSTH = gsMain[indcell, indr+1].subgridspec(2, 1, hspace=0.1, height_ratios=[0.5, 0.5])
            gsRasterPSTH = gsExamples[0, indr].subgridspec(2, 1, hspace=0.1)
            # -- Plot Raster --
            axRaster = plt.subplot(gsRasterPSTH[0])
            pRaster, hcond, zline = extraplots.raster_plot(stfeo, ilet, timeRange,
                                                           trialsEachCond=tec,
                                                           colorEachCond=colorEachCond,
                                                           labels=rasterLabels,
                                                           rasterized=True)
            plt.setp(pRaster, ms=rasterMarkerSize)
            if indr == 1:
                axRaster.set_yticklabels('')
            axRaster.tick_params(labelbottom=False)
            if indcell == 0:
                plt.title(f'{reagentLabels[indr]}', fontweight='normal', fontsize=14)

            # -- Plot PSTH --
            axPSTH = plt.subplot(gsRasterPSTH[1], sharex=axRaster)
            spkCount = spikesanalysis.spiketimes_to_spikecounts(stfeo, ilet, timeVec)
            pPSTH = extraplots.plot_psth(spkCount/binWidth, smoothWinSizePsth, timeVec,
                                         tec, colorEachCond, linestyle=None,
                                         linewidth=lwPsth, downsamplefactor=downsampleFactorPsth,
                                         hidesamples=2)
            plt.xlabel('Time (s)')
            if indr == 0:
                plt.ylabel('Firing rate (spk/s)')
            extraplots.boxoff(axPSTH)
            plt.legend(pPSTH[::-1],['Oddball', 'Standard'], loc='upper left', handlelength=1.5)
            #PSTHyLims = plt.ylim()
            PSTHyLims = [-yMaxEachCell[indcell]/30, yMaxEachCell[indcell]]
            axPSTH.set_ylim(PSTHyLims)
            #axPSTH.set_xlim(timeRangeToShow)
            cstim = studyparams.ODDBALL_SESSION_TYPES[oddballToLoad]['oddball']
            oeiThisCell = dbRow[reagent+cstim+'OEI']
            plt.text(0.2, PSTHyLims[-1]*0.65, f'OEI={oeiThisCell:0.2f}', fontsize=fontSizeLabels)

            stimDuration = bdataList[0]['stimDuration'][-1]
            plot_stim(PSTHyLims, stimDuration)
            plt.show()

if 1:
    # -- Plot summaries --
    OEIsalineEachStim = []
    OEIdoiEachStim = []
    oddballTitles = ['High freq chord', 'Low freq chord', 'Down FM sweep', 'Up FM sweep']
    responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=studyparams.FR_THRESHOLD)
    for inds, stim in enumerate(studyparams.STIMULI):
        if inds in [0,2]:
            gsExamples = gsMain[inds//2, 2].subgridspec(1, 2, wspace=0.15)
        #for inds, stim in enumerate(['high']):
        # -- Select only cells with steady OEI --
        #changeInOEI = (oeindex[:,1]/oeindex[:,0])  # Change from pre to saline
        changeInOEI = celldb['saline'+stim.capitalize()+'OEI'] / celldb['pre'+stim.capitalize()+'OEI']
        steadyOEI = ( (changeInOEI < studyparams.MAX_CHANGE_FACTOR) &
                      (changeInOEI > 1/studyparams.MAX_CHANGE_FACTOR) )
        celldb[stim.capitalize()+'SteadyOEI'] = steadyOEI
        nTrialThreshold = 5
        enoughTrials = ((celldb['saline'+stim.capitalize()+'OddballNtrials']>nTrialThreshold) &
                        (celldb['doi'+stim.capitalize()+'OddballNtrials']>nTrialThreshold))

        #OEIsaline = celldb['saline'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI]
        #OEIdoi = celldb['doi'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI]
        OEIsaline = celldb['saline'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]
        OEIdoi = celldb['doi'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]
        nCellsThisStim = len(OEIsaline)
        medianOEIsaline = np.median(OEIsaline)
        medianOEIdoi = np.median(OEIdoi)

        wstatSvsO, pValSvsO = stats.wilcoxon(OEIsaline, OEIdoi)
        print(f'{stim}\t N={nCellsThisStim} cells\t Median OEI: saline={medianOEIsaline:0.4f}, ' +
              f'DOI={medianOEIdoi:0.4f}   p={pValSvsO:0.4f}')
        print(f'\t DOI>saline: {np.mean(OEIdoi>OEIsaline):0.1%}'+
              f'\t DOI<saline: {np.mean(OEIdoi<OEIsaline):0.1%}')

        #plt.subplot(gsExamples[panelEachStimType[inds][0], panelEachStimType[inds][1]])
        thisAx = plt.subplot(gsExamples[inds%2])
        plt.plot(OEIsaline, OEIdoi, 'o', mec='0.75', mfc='none')
        plt.plot([-1, 1], [-1, 1], 'k-', lw=0.5)
        plt.plot([0, 0], [-1, 1], 'k-', lw=0.5)
        plt.plot([-1, 1], [0, 0], 'k-', lw=0.5)
        plt.plot(medianOEIsaline, medianOEIdoi, '+', color='m', ms=10, mew=2)
        plt.gca().set_aspect('equal', 'box')
        plt.xlim([-1, 1])
        plt.ylim([-1, 1])
        plt.xticks([-1, 0, 1], fontsize=fontSizeTicks)
        plt.yticks([-1, 0, 1], fontsize=fontSizeTicks)
        plt.xlabel('OEI saline', fontsize=fontSizeLabels)
        if inds in[0,2]:
            plt.ylabel('OEI DOI', fontsize=fontSizeLabels)
        else:
            plt.gca().tick_params(labelleft=False)
        plt.title(oddballTitles[inds], fontsize=fontSizeLabels)
        plt.text(0.5, 0.075, f'p = {pValSvsO:0.3f}',
                 transform=thisAx.transAxes, ha='center', fontsize=fontSizeTicks)

        plt.show()

        OEIsalineEachStim.append(OEIsaline)
        OEIdoiEachStim.append(OEIdoi)

        if 0:
            plt.subplot(gsMain[inds//2, inds%2])
            # Plot a bar graph of the OEI for each reagent
            xPos = np.arange(len(studyparams.REAGENTS))
            plt.boxplot([OEIsaline, OEIdoi])
            plt.show()

plt.subplot(gsMain[0,0])
meanOIEsalineEachStim = [np.mean(x) for x in OEIsalineEachStim]
meanOIEDoiEachStim = [np.mean(x) for x in OEIdoiEachStim]
stdOIEsalineEachStim = [np.std(x) for x in OEIsalineEachStim]
stdOIEDoiEachStim = [np.std(x) for x in OEIdoiEachStim]
semOIEsalineEachStim = [stats.sem(x) for x in OEIsalineEachStim]
semOIEDoiEachStim = [stats.sem(x) for x in OEIdoiEachStim]
nCellsEachStim = [len(x) for x in OEIsalineEachStim]
xPos = np.arange(len(studyparams.STIMULI))
xOffset = 0.2
barWidth = 0.25
lineWidth = 2
hbarSal = plt.bar(xPos-xOffset, meanOIEsalineEachStim, width=barWidth, lw=lineWidth,
                  ec=colorsRasterDark['saline'], fc=colorsRasterLight['saline'])
hbarDOI = plt.bar(xPos+xOffset, meanOIEDoiEachStim, width=barWidth, lw=lineWidth,
                  ec=colorsRasterDark['doi'], fc=colorsRasterLight['doi'])
plt.errorbar(xPos-xOffset, meanOIEsalineEachStim, yerr=semOIEsalineEachStim, fmt='none', color='0.5')
plt.errorbar(xPos+xOffset, meanOIEDoiEachStim, yerr=semOIEDoiEachStim, fmt='none', color='0.5')
plt.ylabel('Oddball Enhancement', fontsize=fontSizeLabels)
plt.ylim([0, 0.3])
plt.yticks([0, 0.1, 0.2, 0.3], fontsize=fontSizeTicks)
plt.xlim([-0.5, 3.5])
plt.xticks(xPos, [f'High freq\nchord\n(N={nCellsEachStim[0]})',
                  f'Low freq\nchord\n(N={nCellsEachStim[1]})',
                  f'Down FM \nsweep\n(N={nCellsEachStim[2]})',
                  f'Up FM\nsweep\n(N={nCellsEachStim[3]})'],
           rotation=0, ha='center', fontsize=fontSizeLabels)
starRange = np.array([-xOffset, xOffset])
extraplots.significance_stars(xPos[0]+starRange, 0.31, 0.01, gapFactor=0.2)
extraplots.significance_stars(xPos[1]+starRange, 0.17, 0.01, gapFactor=0.2)
extraplots.significance_stars(xPos[2]+starRange, 0.17, 0.01, gapFactor=0.2)
extraplots.significance_stars(xPos[3]+starRange, 0.17, 0.01, gapFactor=0.2)
extraplots.boxoff(plt.gca())
plt.legend([hbarSal[0], hbarDOI[0]], ['Saline', 'DOI'], loc='upper right', handlelength=1.5, fontsize=fontSizeLabels)
plt.show()

        
'''
gsTopBottom = []
for indsc, stimClass in enumerate(stimClasses):
    gsTopBottom.append(gsMain[indsc].subgridspec(1, 4, hspace=0.1, height_ratios=[0.2, 0.8])
'''

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
