"""
This script plots the evoked responses of each cells under each oddball/standard
condition and under each reagent.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

SAVE_FIGURE = 0
outputDir = '/tmp/'
figFilename = 'plots_oddball_mod_stnd_vs_odd' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.02, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.71]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorOddball = figparams.colors['oddball']
colorStandard = figparams.colors['standard']

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Load celldb of tuning to get steady cells --
maxChangeFactor=1.3
'''
dbTuningFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldbTuning = celldatabase.load_hdf(dbTuningFilename)
#steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steadyParams = ['ToneBaselineFiringRate'] 
steady = studyutils.find_steady_cells(celldbTuning, steadyParams, maxChangeFactor)
'''

# -- Process data --
responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=5)

def oddball_enhancement(evokedStandard, evokedOddball):
    oei = (evokedOddball-evokedStandard) / (evokedOddball+evokedStandard)
    return oei

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Stim types --
stimLabels = ['standard', 'oddball']

# -- Main gridspec --
gsMain = gridspec.GridSpec(len(stimLabels)+1, len(studyparams.STIMULI))
gsMain.update(left=0.07, right=0.98, top=0.92, bottom=0.08, wspace=0.4, hspace=0.3)
axs = []

pColor = '0.5'


for inds, stim in enumerate(studyparams.STIMULI):
    #selectedCells = responsive[stim] & steady
    selectedCells = responsive[stim]
    nCellsThisStim = selectedCells.sum()
    oeindex = np.empty((nCellsThisStim, len(studyparams.REAGENTS)))
    evokedStandard = np.empty((nCellsThisStim, len(studyparams.REAGENTS)))
    evokedOddball = np.empty((nCellsThisStim, len(studyparams.REAGENTS)))
    cstim = stim.capitalize()
    
    for indr, reagent in enumerate(studyparams.REAGENTS):
        evokedStandard[:,indr] = celldb[reagent+cstim+'StandardEvokedFiringRate'][selectedCells]
        evokedOddball[:,indr] = celldb[reagent+cstim+'OddballEvokedFiringRate'][selectedCells]
        
        oeindex[:,indr] = oddball_enhancement(celldb[reagent+cstim+'StandardEvokedFiringRate'][selectedCells],
                                              celldb[reagent+cstim+'OddballEvokedFiringRate'][selectedCells])

    # -- Select only cells with steady OEI --
    if 1:
        changeInOEI = (oeindex[:,1]/oeindex[:,0])  # Change from pre to saline
        steadyOEI = ( (changeInOEI < maxChangeFactor) & (changeInOEI > 1/maxChangeFactor) )
        oeindex = oeindex[steadyOEI, :]
        evokedStandard = evokedStandard[steadyOEI, :]
        evokedOddball = evokedOddball[steadyOEI, :]
        nCellsThisStim = oeindex.shape[0]
    evokedResp = [evokedStandard, evokedOddball]

    modIndexEachLabel = np.empty((nCellsThisStim, len(stimLabels)))
    for indlab, stimLabel in enumerate(stimLabels):

        modIndex1 = studyutils.modulation_index(evokedResp[indlab][:,1], evokedResp[indlab][:,0])
        modIndex2 = studyutils.modulation_index(evokedResp[indlab][:,2], evokedResp[indlab][:,1])
        wstatMI, pValModInd = stats.wilcoxon(modIndex1, modIndex2)
        medianMI1 = np.nanmedian(modIndex1)
        medianMI2 = np.nanmedian(modIndex2)
        print(f'MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')

        modIndexEachLabel[:, indlab] = modIndex2
        
        axs.append(plt.subplot(gsMain[indlab, inds]))
        plt.axhline(0, ls='--', color='0.5', lw=1)
        plt.plot([0, 1, 2], evokedResp[indlab].T, '-', color='0.75')
        plt.plot([0, 1, 2], evokedResp[indlab].T, 'o', mec=pColor, color=pColor)
        for indr, reagent in enumerate(studyparams.REAGENTS):
            medianFR = np.nanmedian(evokedResp[indlab][:,indr])
            plt.plot(indr, medianFR, '_', ms=40, color='0.25')
        plt.xlim([-0.5, 2.5])
        plt.ylabel(f'{stimLabel}', fontsize=fontSizeLabels)
        plt.xticks([0, 1, 2], studyparams.REAGENTS, fontsize=fontSizeTicks)
        wstat1, pVal1 = stats.wilcoxon(evokedResp[indlab][:,0], evokedResp[indlab][:,1])
        wstat2, pVal2 = stats.wilcoxon(evokedResp[indlab][:,1], evokedResp[indlab][:,2])
        plt.text(0.33, 0.95, f'p = {pVal1:0.3f}', transform=axs[-1].transAxes, ha='center',
                 fontsize=fontSizeLabels)
        plt.text(0.66, 0.95, f'p = {pVal2:0.3f}', transform=axs[-1].transAxes, ha='center',
                 fontsize=fontSizeLabels)
        plt.title(f'{cstim} (N={nCellsThisStim}) (p={pValModInd:0.3f})', fontsize=fontSizeLabels)
        plt.show()

    # -- Plot modulation of oddball vs modulation of standard --
    wstatComp, pValStndOdd = stats.wilcoxon(modIndexEachLabel[:,0], modIndexEachLabel[:,1])
    axs.append(plt.subplot(gsMain[2, inds]))
    plt.plot(modIndexEachLabel[:,0], modIndexEachLabel[:,1], 'o', mec='0.75', mfc='none')
    medianMI = np.median(modIndexEachLabel, axis=0)
    plt.plot(medianMI[0], medianMI[1], '+', color='r', ms=10, mew=2)
    plt.xlabel('Mod index standard')
    plt.ylabel('Mod index oddball')
    plt.plot([-1, 1], [-1, 1], 'k-', lw=0.5)
    plt.plot([0, 0], [-1, 1], 'k-', lw=0.5)
    plt.plot([-1, 1], [0, 0], 'k-', lw=0.5)
    plt.gca().set_aspect('equal', 'box')
    plt.xlim([-1, 1])
    plt.ylim([-1, 1])
    plt.title(f'p = {pValStndOdd:0.3f}', fontsize=fontSizeLabels)
    plt.show()

    print('----------')
    print(f'Median change Stnd: {medianMI[0]:0.3f}, Odd: {medianMI[1]:0.3f}  (p = {pValStndOdd:0.3f})')
    print('----------')
sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
