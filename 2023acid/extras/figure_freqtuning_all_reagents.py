"""
This script plots comparisons of frequency tuning under each reagent.
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
figFilename = 'plots_freqtuning_summary' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

if len(sys.argv)==2:
    trialSubset = sys.argv[1]
else:
    trialSubset = ''
if trialSubset not in ['', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.02, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.71]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
#colorOddball = figparams.colors['oddball']
#colorStandard = figparams.colors['standard']

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldbAll = celldatabase.load_hdf(dbFilename)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
    celldb = celldatabase.load_hdf(dbFilename)
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
    celldb = celldatabase.load_hdf(dbFilename)
else:
    celldb = celldbAll
'''
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
if trialSubset == 'running':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_running.h5')
elif trialSubset == 'notrunning':
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning_notrunning.h5')
else:
    dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldb = celldatabase.load_hdf(dbFilename)
'''

# -- Process data --
#maxChangeFactor = 1.2
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR

#responsive = studyutils.find_tone_responsive_cells(celldb, frThreshold=5)
selective = studyutils.find_freq_selective(celldb, minR2=studyparams.MIN_R_SQUARED)
goodFit = studyutils.find_good_gaussian_fit(celldb, minR2=studyparams.MIN_R_SQUARED)
anyFitDOI = ~np.isnan(celldb['doiToneGaussianA'])
posResponse = celldbAll['preToneGaussianA'] > 0
for indr, reagent in enumerate(studyparams.REAGENTS):
    celldb[reagent+'ToneGaussianMax'] = ( celldb[reagent+'ToneGaussianA'] +
                                          celldb[reagent+'ToneGaussianY0'] )
#steadyParams = ['ToneBaselineFiringRate', 'ToneGaussianSigma', 'ToneGaussianMax'] 
#steadyParams = ['ToneBaselineFiringRate', 'ToneGaussianMax'] 
#steadyParams = ['ToneBaselineFiringRate', 'ToneGaussianSigma'] 
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldb, steadyParams, maxChangeFactor)

'''
steady = studyutils.find_tone_steady_cells(celldb, maxChangeFactor=maxChangeFactor)
selective = studyutils.find_freq_selective(celldb, minR2=studyparams.MIN_R_SQUARED)
posResponse = celldb['preToneGaussianA'] > 0
#steadyParam = studyutils.find_steady_by_param(celldb, 'ToneGaussianSigma', maxChangeFactor=1.2)
steadyA = studyutils.find_steady_by_param(celldb, 'ToneGaussianA', maxChangeFactor=1.2)
steadyBaseline = studyutils.find_steady_by_param(celldb, 'ToneBaselineFiringRate', maxChangeFactor=1.2)
steadySigma = studyutils.find_steady_by_param(celldb, 'ToneGaussianSigma', maxChangeFactor=1.2)

steadyA = studyutils.find_steady_by_param(celldb, 'ToneGaussianA', maxChangeFactor=1.2)
steadyBaseline = studyutils.find_steady_cells(celldb, 'ToneBaselineFiringRate',
                                              maxChangeFactor=maxChangeFactor)
steadySigma = studyutils.find_steady_cells(celldb, 'ToneGaussianSigma',
                                           maxChangeFactor=maxChangeFactor)
#steadyGaussianMax = studyutils.find_steady_cells(celldb, 'ToneGaussianMax',
#                                                 maxChangeFactor=maxChangeFactor)

'''



#metrics = ['ToneBaselineFiringRate', 'ToneFiringRateBestFreq', 'ToneAvgEvokedFiringRate']
metrics = ['ToneGaussianA', 'ToneGaussianX0', 'ToneGaussianSigma', 'ToneGaussianY0',
           'ToneGaussianMax']

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

# -- Main gridspec --
#gsMain = gridspec.GridSpec(1, len(studyparams.REAGENTS))
gsMain = gridspec.GridSpec(1, len(metrics))
gsMain.update(left=0.07, right=0.98, top=0.92, bottom=0.08, wspace=0.4, hspace=0.3)
axs = []

pColor = '0.5'

for indm, metric in enumerate(metrics):

    thisMetricFiring = np.empty((len(celldb), len(studyparams.REAGENTS)))
    for indr, reagent in enumerate(studyparams.REAGENTS):
        #baselineFiring[:, indr] = celldb[reagent+'ToneBaselineFiringRate']
        #baselineFiring[:, indr] = celldb[reagent+'ToneAvgEvokedFiringRate']
        #if metric=='maxPredFiringRate':
        #    thisMetricFiring[:, indr] = celldb[reagent+'ToneGaussianA']+celldb[reagent+'ToneGaussianY0']
        #else:
        thisMetricFiring[:, indr] = celldb[reagent + metric]

    #thisMetricFiring = thisMetricFiring[responsive & steady & selective, :] 
    #thisMetricFiring = thisMetricFiring[responsive & steady & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[responsive & steady & steadyParam & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[responsive & steady & steadyBaseline & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[responsive & selective & steadyGaussianMax & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[responsive & steadyA & steadySigma & steadyBaseline & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[responsive & posResponse & steady, :] 
    #thisMetricFiring = thisMetricFiring[selective & posResponse & steady, :] 
    #thisMetricFiring = thisMetricFiring[goodFit & posResponse, :] 
    #thisMetricFiring = thisMetricFiring[goodFit & posResponse & steady & anyFitDOI, :] 
    #thisMetricFiring = thisMetricFiring[responsive & goodFit & posResponse & steady, :] 
    #thisMetricFiring = thisMetricFiring[goodFit & posResponse & anyFitDOI, :] 
    thisMetricFiring = thisMetricFiring[goodFit & posResponse & steady & anyFitDOI, :] 
    
    modIndex1 = studyutils.modulation_index(thisMetricFiring[:,1], thisMetricFiring[:,0])
    modIndex2 = studyutils.modulation_index(thisMetricFiring[:,2], thisMetricFiring[:,1])
    wstatMI, pValModInd = stats.wilcoxon(modIndex1, modIndex2)
    medianMI1 = np.nanmedian(modIndex1)
    medianMI2 = np.nanmedian(modIndex2)
    print(f'MI: {medianMI1:0.3f}, {medianMI2:0.3f} (p = {pValModInd:0.3f})')
    
    axs.append(plt.subplot(gsMain[indm]))
    plt.axhline(0, ls='--', color='0.5', lw=1)
    plt.plot([0, 1, 2], thisMetricFiring.T, '-', color='0.75')
    plt.plot([0, 1, 2], thisMetricFiring.T, 'o', mec=pColor, color=pColor)
    for indr, reagent in enumerate(studyparams.REAGENTS):
        medianFR = np.nanmedian(thisMetricFiring[:,indr])
        plt.plot(indr, medianFR, '_', ms=40, color='0.25')
    plt.xlim([-0.5, 2.5])
    plt.ylabel(f'{metric} (spk/s)', fontsize=fontSizeLabels)
    plt.xticks([0, 1, 2], studyparams.REAGENTS, fontsize=fontSizeTicks)
    wstat1, pVal1 = stats.wilcoxon(thisMetricFiring[:,0], thisMetricFiring[:,1])
    wstat2, pVal2 = stats.wilcoxon(thisMetricFiring[:,1], thisMetricFiring[:,2])
    plt.text(0.33, 0.95, f'p = {pVal1:0.3f}', transform=axs[-1].transAxes, ha='center',
             fontsize=fontSizeLabels)
    plt.text(0.66, 0.95, f'p = {pVal2:0.3f}', transform=axs[-1].transAxes, ha='center',
             fontsize=fontSizeLabels)
    #plt.title(f'{cstim} (N={nCellsThisStim})', fontsize=fontSizeLabels)
    plt.title(f'N = {len(thisMetricFiring)} (p={pValModInd:0.3f})', fontsize=fontSizeLabels)
    plt.show()

sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
