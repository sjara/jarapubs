"""
Show the change in oddball compared to standard when changing from saline to DOI.
"""

import os
import sys
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


stimDuration = 0.1 ### Only used for the plot (not analysis), BUT THIS SHOULD NOT BE HARDCODED!

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'figure_modulation_odd_vs_stnd' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [14, 3] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

#labelPosX = [0.01, 0.30, 0.49, 0.66, 0.83]   # Horiz position for panel labels
labelPosX = [0.005, 0.30, 0.66]   # Horiz position for panel labels
labelPosY = [0.92]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorsPSTH = [[figparams.colorsLight['saline'], figparams.colorsLight['doi']],
              [figparams.colors['saline'], figparams.colors['doi']]]
pColor = '0.5'

if 1:
    cellToPlot = 46  #[86,  46]  #  High:84, 86  Low: 1065
    oddballToPlot = 'FM_Up'  #['HighFreq', 'FM_Up']
    yMaxEachCell = 200  #[80, 180]
else:
    cellToPlot = 86  #[86,  46]  #  High:84, 86  Low: 1065
    oddballToPlot = 'HighFreq'  #, 'FM_Up']
    yMaxEachCell = 80  #[80, 180]
'''
    cellToPlot = 86  #[86,  46]  #  High:84, 86  Low: 1065
    oddballToPlot = 'HighFreq'  #, 'FM_Up']
    yMaxEachCell = 80  #[80, 180]
'''

# -- Raster and PSTH parameters --
timeRange = [-0.2, 0.35]
binWidth = 0.010
timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
smoothWinSizePsth = 2 
lwPsth = 2.5
downsampleFactorPsth = 1

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Calculate modulation index Saline-to-DOI for standards and for oddballs --
for inds, stim in enumerate(studyparams.STIMULI):
    cstim = stim.capitalize()
    evokedStandard = np.empty((len(celldb), len(studyparams.REAGENTS)))
    evokedOddball = np.empty((len(celldb), len(studyparams.REAGENTS)))
    for indr, reagent in enumerate(studyparams.REAGENTS):
        evokedStandard[:,indr] = celldb[reagent+cstim+'StandardEvokedFiringRate']
        evokedOddball[:,indr] = celldb[reagent+cstim+'OddballEvokedFiringRate']
        celldb[reagent+cstim+'OEI'] = studyutils.modulation_index(evokedOddball[:,indr], evokedStandard[:,indr])
    # -- Calculate modulation index between DOI (ind=2) and Saline (ind=1) --
    modIndexStandard = studyutils.modulation_index(evokedStandard[:,2], evokedStandard[:,1])
    modIndexOddball = studyutils.modulation_index(evokedOddball[:,2], evokedOddball[:,1])
    celldb[cstim+'ModIndexStandard'] = modIndexStandard
    celldb[cstim+'ModIndexOddball'] = modIndexOddball
    

def plot_stim(yLims, stimDuration, stimLineWidth=6, stimColor=figparams.colorStim):
    yPos = 0.85*yLims[-1] + 0.075*(yLims[-1]-yLims[0])
    pstim = plt.plot([0, stimDuration], 2*[yPos], lw=stimLineWidth, color=stimColor,
                     clip_on=False, solid_capstyle='butt')
    return pstim[0]

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

stimClasses = ['FM', 'pureTone']


# -- Main gridspec --
#gsMain = gridspec.GridSpec(1, 5, width_ratios=[0.4, 0.15, 0.15, 0.15, 0.15])
gsMain = gridspec.GridSpec(1, 5, width_ratios=[0.28]+[0.18]*4)
gsMain.update(left=0.09, right=0.98, top=0.92, bottom=0.15, wspace=0.2, hspace=0.4)
gsRasterPSTHs = gsMain[0].subgridspec(2, 1, hspace=0.1)

# -- Plot PSTHs --
dbRow = celldb.iloc[cellToPlot]
oneCell = ephyscore.Cell(dbRow)
reagentsToPlot = ['saline', 'doi']
reagentLabels = ['Saline', 'DOI']
oddballToLoad = oddballToPlot
cstim = studyparams.ODDBALL_SESSION_TYPES[oddballToLoad]['oddball']
modIndexEachStimCond = [dbRow[cstim+'ModIndexStandard'],
                        dbRow[cstim+'ModIndexOddball']]
                        
axsPSTH = []
ephysEachReagent = {}
ephysEachStimCond = {}
for indr, reagent in enumerate(reagentsToPlot):
    (stfeo, ilet, tec, clabels, bdataList) = \
        studyutils.load_oddball_data_one_cell(oneCell, oddballToLoad,
                                              reagent, timeRange, nPreOdd=8)
    ephysEachReagent[reagent] = [stfeo, ilet, tec, clabels]

(stfeo, ilet, tec) = studyutils.combine_index_limits(ephysEachReagent['saline'][0],
                                                     ephysEachReagent['doi'][0],
                                                     ephysEachReagent['saline'][1],
                                                     ephysEachReagent['doi'][1])
trialsEachCond = studyutils.combine_trials_each_cond(ephysEachReagent['saline'][2],
                                                     ephysEachReagent['doi'][2])
trialsEachCondLabels = ['saline_standard', 'saline_oddball',
                        'doi_standard', 'doi_oddball']
spkCount = spikesanalysis.spiketimes_to_spikecounts(stfeo, ilet, timeVec)
#clabels
#stimToPlot = ['standard', 'oddball']
                        
for indsc, stimcond in enumerate(clabels):
    # -- Plot PSTH --
    axPSTH = plt.subplot(gsRasterPSTHs[indsc])
    extraplots.adjust_pos(axPSTH, [-0.03, 0, 0, 0])
    colorEachCond = colorsPSTH[indsc]
    condsToCompare = np.array([0, 2]) + indsc
    trialsTheseConds= trialsEachCond[:,condsToCompare]
    pPSTH = extraplots.plot_psth(spkCount/binWidth, smoothWinSizePsth, timeVec,
                                 trialsTheseConds, colorEachCond, linestyle=None,
                                 linewidth=lwPsth, downsamplefactor=downsampleFactorPsth)
    plt.ylabel('Firing rate\n(spk/s)', fontsize=fontSizeLabels)
    extraplots.boxoff(axPSTH)
    plt.legend(pPSTH, reagentLabels, loc='upper left', handlelength=1.5)
    #PSTHyLims = plt.ylim()
    PSTHyLims = [-yMaxEachCell/30, yMaxEachCell]
    axPSTH.set_ylim(PSTHyLims)
    axPSTH.set_xlim(timeRange)
    if indsc==0:
        axPSTH.tick_params(labelbottom=False)
    else:
        plt.xlabel('Time (s)', fontsize=fontSizeLabels)
    plt.text(0.175, PSTHyLims[-1]*0.75, f'{stimcond.capitalize()}', fontweight='bold', fontsize=fontSizeLabels)
    plt.text(0.175, PSTHyLims[-1]*0.55, f'MI = {modIndexEachStimCond[indsc]:0.2f}', fontsize=fontSizeLabels)

    plot_stim(PSTHyLims, stimDuration)

# -- Show panel labels --
#for indplabel, plabel in enumerate(['A','B','C','D','E']):
for indplabel, plabel in enumerate(['A','B','C']):
    plt.gca().annotate(plabel, xy=(labelPosX[indplabel],labelPosY[0]), xycoords='figure fraction',
                 fontsize=fontSizePanel, fontweight='bold')


# -- Summary for each stim type --
oddballTitles = ['High freq chord', 'Low freq chord', 'Down FM sweep', 'Up FM sweep']
responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=studyparams.FR_THRESHOLD)
for inds, stim in enumerate(studyparams.STIMULI):
    changeInOEI = celldb['saline'+stim.capitalize()+'OEI'] / celldb['pre'+stim.capitalize()+'OEI']
    steadyOEI = ( (changeInOEI < studyparams.MAX_CHANGE_FACTOR) &
                  (changeInOEI > 1/studyparams.MAX_CHANGE_FACTOR) )
    celldb[stim.capitalize()+'SteadyOEI'] = steadyOEI

    #[responsive[stim] & steadyOEI]
    modIndexStandard = celldb[stim.capitalize()+'ModIndexStandard'][responsive[stim] & steadyOEI]
    modIndexOddball = celldb[stim.capitalize()+'ModIndexOddball'][responsive[stim] & steadyOEI]
    medianMIstnd = np.median(modIndexStandard)
    medianMIodd = np.median(modIndexOddball)

    thisAx = plt.subplot(gsMain[0, inds+1])
    #print(thisAx.get_position())
    
    wstatComp, pValStndOdd = stats.wilcoxon(modIndexStandard, modIndexOddball)
    plt.plot(modIndexStandard, modIndexOddball, 'o', mec='0.75', mfc='none')
    plt.plot(medianMIstnd, medianMIodd, '+', color='m', ms=10, mew=2)
    plt.xlabel('MI standard', fontsize=fontSizeLabels)
    if inds in [0,2]:
        plt.ylabel('MI oddball', fontsize=fontSizeLabels)
    else:
        thisAx.tick_params(labelleft=False)
        extraplots.adjust_pos(thisAx, [-0.01, 0, 0, 0])
    if inds in [2,3]:
        extraplots.adjust_pos(thisAx, [0.02, 0, 0, 0])

    plt.plot([-1, 1], [-1, 1], 'k-', lw=0.5)
    plt.plot([0, 0], [-1, 1], 'k-', lw=0.5)
    plt.plot([-1, 1], [0, 0], 'k-', lw=0.5)
    plt.gca().set_aspect('equal', 'box')
    plt.xlim([-1, 1])
    plt.ylim([-1, 1])
    plt.xticks([-1, 0, 1], fontsize=fontSizeTicks)
    plt.yticks([-1, 0, 1], fontsize=fontSizeTicks)
    #plt.title(f'p = {pValStndOdd:0.6f}', fontsize=fontSizeLabels)
    plt.title(oddballTitles[inds], fontsize=fontSizeLabels)
    plt.text(0.5, 0.9, f'p = {pValStndOdd:0.3f}',
             transform=thisAx.transAxes, ha='center', fontsize=fontSizeTicks)
    plt.show()

    print(f'{stim}: Median change Stnd: {medianMIstnd:0.3f}, Odd: {medianMIodd:0.3f}  (p = {pValStndOdd:0.6f})')
    print(f'\t oddball>standard: {np.mean(modIndexOddball>modIndexStandard):0.1%}'+
        f'\t oddball<standard: {np.mean(modIndexOddball<modIndexStandard):0.1%}')

    plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
