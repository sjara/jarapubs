"""
Show the oddball enhancement index (OEI) for either early or late days.
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



SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'figure_early_vs_late' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [14, 4] # In inches

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

#labelPosX = [0.005, 0.34, 0.68]   # Horiz position for panel labels
labelPosX = [0.005, 0.27, 0.515, 0.76]   # Horiz position for panel labels
labelPosY = [0.95]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorsRasterDark = figparams.colors
colorsRasterLight = figparams.colorsLight
pColor = '0.5'
rasterMarkerSize = 0.5


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

# -- Select only early or late days --
earlyDatesInds = [0,1]  # First two days
lateDatesInds = [2,3,4]   # Last two days
cellsEarly = np.zeros(len(celldb), dtype=bool)
cellsLate = np.zeros(len(celldb), dtype=bool)
for subject in celldb.subject.unique():
    subjectCells = celldb[celldb.subject==subject]
    datesThisSubject = subjectCells.date.unique()
    print(f'{subject}: {datesThisSubject}')
    for date in datesThisSubject:
        if date in datesThisSubject[earlyDatesInds]:
            cellsEarly = cellsEarly | ((celldb.subject==subject) & (celldb.date==date))
        elif date in datesThisSubject[lateDatesInds]:
            cellsLate = cellsLate | ((celldb.subject==subject) & (celldb.date==date))
        else:
            pass
celldbEarly = celldb.loc[cellsEarly].copy()
celldbLate = celldb.loc[cellsLate].copy()

daysLabels = ['All days', 'First two days', 'Last three days', 'Last three days (saline matched)']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

stimClasses = ['FM', 'pureTone']
# -- Main gridspec --
gsMain = gridspec.GridSpec(1, len(daysLabels))
gsMain.update(left=0.04, right=0.99, top=0.94, bottom=0.16, wspace=0.15)

# -- Show panel labels --
for indp, plabel in enumerate(['A','B','C', 'D']):
    plt.figtext(labelPosX[indp], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')

pValoffsety = [0.05, 0.04, 0.05, 0.035]

OEIsalineEachStimEachPeriod = []
OEIdoiEachStimEachPeriod = []

# -- Loops through ranges of days --
for indperiod, celldb in enumerate([celldbAll, celldbEarly, celldbLate, celldbLate]):

    # -- Calculate oddball enhancement index (OEI) --
    for inds, stim in enumerate(studyparams.STIMULI):
        cstim = stim.capitalize()
        oeindex = np.empty((len(celldb), len(studyparams.REAGENTS)))
        for indr, reagent in enumerate(studyparams.REAGENTS):
            standardEvokedFR = celldb[reagent+cstim+'StandardEvokedFiringRate']
            oddballEvokedFR = celldb[reagent+cstim+'OddballEvokedFiringRate']
            celldb[reagent+cstim+'OEI'] = studyutils.modulation_index(oddballEvokedFR, standardEvokedFR)
            #oeindex[:,indr] = studyutils.modulation_index(oddballEvokedFR, standardEvokedFR)

    OEIsalineEachStim = []
    OEIdoiEachStim = []
    pValSvsOperStim = []
    nCellsEachStimEachPeriod = []
    nUniqueCellsEachStimEachPeriod = []
    oddballTitles = ['High freq chord', 'Low freq chord', 'Down FM sweep', 'Up FM sweep']
    responsive = studyutils.find_oddball_responsive_cells(celldb, frThreshold=studyparams.FR_THRESHOLD)
    for inds, stim in enumerate(studyparams.STIMULI):
        changeInOEI = celldb['saline'+stim.capitalize()+'OEI'] / celldb['pre'+stim.capitalize()+'OEI']
        steadyOEI = ( (changeInOEI < studyparams.MAX_CHANGE_FACTOR) &
                      (changeInOEI > 1/studyparams.MAX_CHANGE_FACTOR) )
        celldb[stim.capitalize()+'SteadyOEI'] = steadyOEI
        nTrialThreshold = 5
        enoughTrials = ((celldb['saline'+stim.capitalize()+'OddballNtrials']>nTrialThreshold) &
                        (celldb['doi'+stim.capitalize()+'OddballNtrials']>nTrialThreshold))

        OEIsaline = celldb['saline'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]
        OEIdoi = celldb['doi'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]

        OEIsaline = np.array(OEIsaline)
        OEIdoi = np.array(OEIdoi)
        
        # -- Match saline OEI for late period to that of early period --
        if indperiod==3:
            np.random.seed(10)
            #OEIsalineEarly = OEIsalineEachStimEachPeriod[1][inds].to_numpy()
            OEIsalineEarly = OEIsalineEachStimEachPeriod[1][inds]
            meanOEIsalineEarly = np.mean(OEIsalineEarly)
            #OEIsaline = np.array(OEIsaline)
            #OEIdoi = np.array(OEIdoi)
            minSize = min(len(OEIsaline), len(OEIsalineEarly))
            #newOEIind = np.random.choice(len(OEIsaline), size=len(OEIsalineEarly), replace=True)
            newOEIind = np.random.choice(len(OEIsaline), size=minSize, replace=False)
            #newOEIind = np.arange(len(OEIsaline))
            newOEI = OEIsaline[newOEIind]

            for ind in range(1000):  # Max iterations to prevent infinite loop
                meanNewOEI = np.mean(newOEI)
                if np.isclose(meanNewOEI, meanOEIsalineEarly, atol=0.005):  # Stop if close enough
                    print(f'iter: {ind}')
                    break
                if meanNewOEI < meanOEIsalineEarly:
                    idxToReplace = np.argmin(newOEI)
                else:
                    idxToReplace = np.argmax(newOEI)
                #newOEIind[idxToReplace] = np.random.choice(range(len(OEIsaline)))
                newOEIind = np.delete(newOEIind, idxToReplace)
                newOEI = OEIsaline[newOEIind]
            OEIsaline = newOEI
            OEIdoi = OEIdoi[newOEIind]
            nCellsThisStim = len(OEIsaline)
            nUniqueCellsThisStim = len(np.unique(newOEIind))
            nCellsEachStimEachPeriod.append(nCellsThisStim)
            nUniqueCellsEachStimEachPeriod.append(nUniqueCellsThisStim)
       
        medianOEIsaline = np.median(OEIsaline)
        medianOEIdoi = np.median(OEIdoi)
        OEIsalineEachStim.append(OEIsaline)
        OEIdoiEachStim.append(OEIdoi)
        wstatSvsO, pValSvsO = stats.wilcoxon(OEIsaline, OEIdoi)
        pValSvsOperStim.append(pValSvsO)

    ax = plt.subplot(gsMain[indperiod])
    #print(ax.get_position())
    print(pValSvsOperStim)

    OEIsalineEachStimEachPeriod.append(OEIsalineEachStim)
    OEIdoiEachStimEachPeriod.append(OEIdoiEachStim)
    
    meanOIEsalineEachStim = [np.mean(x) for x in OEIsalineEachStim]
    meanOIEDoiEachStim = [np.mean(x) for x in OEIdoiEachStim]

    # -- Use median instead --
    if 0:
        meanOIEsalineEachStim = [np.median(x) for x in OEIsalineEachStim]
        meanOIEDoiEachStim = [np.median(x) for x in OEIdoiEachStim]
    
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
    plt.axhline(0, color='k', lw=1)
    plt.errorbar(xPos-xOffset, meanOIEsalineEachStim, yerr=semOIEsalineEachStim, fmt='none', color='0.5')
    plt.errorbar(xPos+xOffset, meanOIEDoiEachStim, yerr=semOIEDoiEachStim, fmt='none', color='0.5')
    if indperiod==0:
        plt.ylabel('Oddball Enhancement', fontsize=fontSizeLabels)
    plt.ylim([-0.02, 0.34])
    plt.yticks([0, 0.1, 0.2, 0.3], fontsize=fontSizeTicks)
    if indperiod!=0:
        plt.gca().tick_params(labelleft=False)
    plt.xlim([-0.5, 3.5])
    '''
    plt.xticks(xPos, [f'High freq\nchord\n(N={nCellsEachStim[0]})',
                      f'Low freq\nchord\n(N={nCellsEachStim[1]})',
                      f'Down FM \nsweep\n(N={nCellsEachStim[2]})',
                      f'Up FM\nsweep\n(N={nCellsEachStim[3]})'],
               rotation=0, ha='center', fontsize=fontSizeLabels)
    '''
    plt.xticks(xPos, [f'High F\nchord\n(N={nCellsEachStim[0]})',
                      f'Low F\nchord\n(N={nCellsEachStim[1]})',
                      f'Down FM \nsweep\n(N={nCellsEachStim[2]})',
                      f'Up FM\nsweep\n(N={nCellsEachStim[3]})'],
               rotation=0, ha='center', fontsize=fontSizeLabels)
    starRange = np.array([-xOffset, xOffset])
    
    '''
    extraplots.significance_stars(xPos[0]+starRange, 0.31, 0.01, gapFactor=0.2)
    extraplots.significance_stars(xPos[1]+starRange, 0.17, 0.01, gapFactor=0.2)
    extraplots.significance_stars(xPos[2]+starRange, 0.17, 0.01, gapFactor=0.2)
    extraplots.significance_stars(xPos[3]+starRange, 0.17, 0.01, gapFactor=0.2)
    '''
    for indst, pVal in enumerate(pValSvsOperStim):
        if indst==0:
            pValypos = meanOIEsalineEachStim[0]+pValoffsety[indperiod]
        else:
            pValypos = 0.23
        fontweight = 'bold' if pVal<0.05 else 'normal'
        plt.text(xPos[indst], pValypos, f'{pVal:0.3f}', fontsize=fontSizeTicks, ha='center',
                 fontweight=fontweight, rotation=0)
    extraplots.boxoff(plt.gca())
    plt.legend([hbarSal[0], hbarDOI[0]], ['Saline', 'DOI'], loc='upper right', handlelength=1.5, fontsize=fontSizeLabels)
    plt.title(daysLabels[indperiod], fontsize=fontSizeLabels)
    plt.show()


# -- Statistics between periods --
print('\nComparison between early and late (not matched)')
latePeriodInd = 2
for inds, stim in enumerate(studyparams.STIMULI):
    earlyEffect = OEIsalineEachStimEachPeriod[1][inds]-OEIdoiEachStimEachPeriod[1][inds]
    lateEffect = OEIsalineEachStimEachPeriod[latePeriodInd][inds]-OEIdoiEachStimEachPeriod[latePeriodInd][inds]
    ustatEvsL, pValEvsL = stats.mannwhitneyu(earlyEffect, lateEffect)
    print(f'{stim}: {pValEvsL:0.4f} \t Median early: {np.median(earlyEffect):0.3f} \t Median late: {np.median(lateEffect):0.3f}')

print('\nComparison between early and late (matched)')
latePeriodInd = 3
for inds, stim in enumerate(studyparams.STIMULI):
    earlyEffect = OEIsalineEachStimEachPeriod[1][inds]-OEIdoiEachStimEachPeriod[1][inds]
    lateEffect = OEIsalineEachStimEachPeriod[latePeriodInd][inds]-OEIdoiEachStimEachPeriod[latePeriodInd][inds]
    ustatEvsL, pValEvsL = stats.mannwhitneyu(earlyEffect, lateEffect)
    print(f'{stim}: {pValEvsL:0.4f} \t Median early: {np.median(earlyEffect):0.3f} \t Median late: {np.median(lateEffect):0.3f}')

print('')
print('nCells:',nCellsEachStimEachPeriod)
print('nUnique:',nUniqueCellsEachStimEachPeriod)

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
