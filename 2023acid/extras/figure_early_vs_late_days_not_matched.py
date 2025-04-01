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
figFilename = 'plots_early_vs_late' # Do not include extension
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

labelPosX = [0.005, 0.34, 0.68]   # Horiz position for panel labels
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

daysLabels = ['All days', 'First two days', 'Last three days']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

stimClasses = ['FM', 'pureTone']
# -- Main gridspec --
gsMain = gridspec.GridSpec(1, len(daysLabels))
gsMain.update(left=0.05, right=0.99, top=0.94, bottom=0.16, wspace=0.3)

# -- Show panel labels --
for indp, plabel in enumerate(['A','B','C']):
    plt.figtext(labelPosX[indp], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')

pValoffsety = [0.05, 0.04, 0.05]

OEIsalineEachStimEachPeriod = []
OEIdoiEachStimEachPeriod = []

# -- Loops through ranges of days --
for indperiod, celldb in enumerate([celldbAll, celldbEarly, celldbLate]):

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

        # -- Select only cells with high OEI --
        if 0:
            #enoughTrials = enoughTrials & (celldb['saline'+stim.capitalize()+'OEI']>0)
            enoughTrials = enoughTrials & (celldb['doi'+stim.capitalize()+'OEI']>0)
        # -- Make saline OEI for late period have values over the min value of early period --
        if 0:
            if indperiod==2:
                minOEI = np.min(OEIsalineEachStimEachPeriod[1][inds])
                enoughTrials = enoughTrials & (celldb['saline'+stim.capitalize()+'OEI']>minOEI)
        
        OEIsaline = celldb['saline'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]
        OEIdoi = celldb['doi'+stim.capitalize()+'OEI'][responsive[stim] & steadyOEI & enoughTrials]

        # -- Match saline OEI for late period to that of early period --
        if 0:
            np.random.seed(10)
            if indperiod==2:
                OEIsalineEarly = OEIsalineEachStimEachPeriod[1][inds].to_numpy()
                meanOEIsalineEarly = np.mean(OEIsalineEarly)
                OEIsaline = np.array(OEIsaline)
                OEIdoi = np.array(OEIdoi)
                minSize = min(len(OEIsaline), len(OEIsalineEarly))
                #newOEIind = np.random.choice(len(OEIsaline), size=len(OEIsalineEarly), replace=True)
                newOEIind = np.random.choice(len(OEIsaline), size=minSize, replace=False)
                newOEI = OEIsaline[newOEIind]

                for ind in range(1000):  # Max iterations to prevent infinite loop
                    meanNewOEI = np.mean(newOEI)
                    if np.isclose(meanNewOEI, meanOEIsalineEarly, atol=0.01):  # Stop if close enough
                        print(f'iter: {ind}')
                        break
                    if meanNewOEI < meanOEIsalineEarly:
                        idxToReplace = np.argmin(newOEI)
                    else:
                        idxToReplace = np.argmax(newOEI)
                    newOEIind[idxToReplace] = np.random.choice(range(len(OEIsaline)))
                    newOEI = OEIsaline[newOEIind]
                OEIsaline = newOEI
                OEIdoi = OEIdoi[newOEIind]
       
        nCellsThisStim = len(OEIsaline)
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
    plt.ylabel('Oddball Enhancement', fontsize=fontSizeLabels)
    plt.ylim([-0.02, 0.34])
    plt.yticks([0, 0.1, 0.2, 0.3], fontsize=fontSizeTicks)
    plt.xlim([-0.5, 3.5])
    plt.xticks(xPos, [f'High freq\nchord\n(N={nCellsEachStim[0]})',
                      f'Low freq\nchord\n(N={nCellsEachStim[1]})',
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
            pValypos = 0.22
        fontweight = 'bold' if pVal<0.05 else 'normal'
        plt.text(xPos[indst], pValypos, f'p={pVal:0.3f}', fontsize=fontSizeTicks, ha='center',
                 fontweight=fontweight, rotation=0)
    extraplots.boxoff(plt.gca())
    plt.legend([hbarSal[0], hbarDOI[0]], ['Saline', 'DOI'], loc='upper right', handlelength=1.5, fontsize=fontSizeLabels)
    plt.title(daysLabels[indperiod], fontsize=fontSizeLabels)
    plt.show()


if 0:
    # -- Main gridspec --
    gsMain = gridspec.GridSpec(len(stimClasses), 3, width_ratios=[0.16, 0.52, 0.32])
    gsMain.update(left=0.07, right=0.99, top=0.94, bottom=0.08, wspace=0.25, hspace=0.4)

    # -- Show panel labels --
    for indp, plabel in enumerate(['A','B','C','D']):
        plt.figtext(labelPosX[indp], labelPosY[0], plabel, fontsize=fontSizePanel, fontweight='bold')
    for indp, plabel in enumerate(['E','F','G','H']):
        plt.figtext(labelPosX[indp], labelPosY[1], plabel, fontsize=fontSizePanel, fontweight='bold')

    panelEachStimType = [[0,3], [0,4], [1,3], [1,4]]

    # -- Plot summaries --
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
    

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
