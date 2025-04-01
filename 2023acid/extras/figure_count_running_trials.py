"""
Show a comparison of how many trials the animal was running vs non-running.
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

pColor = '0.5'
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
dbFilenameRun = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
dbFilenameNotrun = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
condLabels = ['notrunning', 'running']
celldbAll = celldatabase.load_hdf(dbFilename)
celldbs = [celldatabase.load_hdf(dbFilenameNotrun),
           celldatabase.load_hdf(dbFilenameRun)]


# -- Estimate how many trials in each session under each running condition. --
nTrialsEachRunningCond = [[], []]
for inds, subject in enumerate(studyparams.SUBJECTS):
    for indc, celldbThisCond in enumerate(celldbs):
        celldbThisSubject = celldbThisCond.query('subject==@subject')
        datesThisSubject = list(celldbThisSubject.date.unique())
        nTrialsThisRunning = []
        for oneDate in datesThisSubject:
            celldbThisDate = celldbThisSubject.query('date==@oneDate')
            nTrialsThisDate = [celldbThisDate.preToneNtrials.iloc[0],
                               celldbThisDate.salineToneNtrials.iloc[0],
                               celldbThisDate.doiToneNtrials.iloc[0]]
            nTrialsEachRunningCond[indc].append(nTrialsThisDate)
            #print(f'{subject} {oneDate} {condLabels[indc]}: {nTrialsThisDate}')

# -- Fraction of running trials under each reagent condition. --
fractionRunning = ( np.array(nTrialsEachRunningCond[1]) /
                    (np.array(nTrialsEachRunningCond[0])+np.array(nTrialsEachRunningCond[1])) )
pcRunning = 100*fractionRunning

modIndex1 = studyutils.modulation_index(pcRunning[:,1], pcRunning[:,0])
modIndex2 = studyutils.modulation_index(pcRunning[:,2], pcRunning[:,1])
wstatMI, pValModInd = stats.wilcoxon(modIndex1, modIndex2)
medianMI1 = np.nanmedian(modIndex1)
medianMI2 = np.nanmedian(modIndex2)
print(f'MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')

plt.clf()
axs = [plt.gca()]
#plt.axhline(0, ls='--', color='0.5', lw=1)
plt.plot([0, 1, 2], pcRunning.T, '-', color='0.75')
plt.plot([0, 1, 2], pcRunning.T, 'o', mec=pColor, color=pColor)
for indr, reagent in enumerate(studyparams.REAGENTS):
    medianFR = np.nanmedian(pcRunning[:,indr])
    plt.plot(indr, medianFR, '_', ms=40, color='0.25')
plt.xlim([-0.5, 2.5])
plt.ylim([0, 100])
plt.ylabel(f'Percent trials running (%)', fontsize=fontSizeLabels)
plt.xticks([0, 1, 2], studyparams.REAGENTS, fontsize=fontSizeTicks)
wstat1, pVal1 = stats.wilcoxon(pcRunning[:,0], pcRunning[:,1])
wstat2, pVal2 = stats.wilcoxon(pcRunning[:,1], pcRunning[:,2])
plt.text(0.33, 0.95, f'p = {pVal1:0.3f}', transform=axs[-1].transAxes, ha='center',
         fontsize=fontSizeLabels)
plt.text(0.66, 0.95, f'p = {pVal2:0.3f}', transform=axs[-1].transAxes, ha='center',
         fontsize=fontSizeLabels)
#plt.title(f'{cstim} (N={nCellsThisStim})', fontsize=fontSizeLabels)
plt.title(f'N = {len(pcRunning)} (p={pValModInd:0.3f})', fontsize=fontSizeLabels)
plt.show()



