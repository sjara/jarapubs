"""
Analyze oddball sessions under each reagent.

This script adds columns to the database with the firing rate during baseline and evoked periods.

NOTE: when estimating the response to the standard, we use the last 4 presentations before each oddball.

In Santiago's workstation it took between 5-7min to run.
In Santiago's laptop it took between 11min to run.
"""

import os
import sys
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
from jaratoolbox import extraplots
from jaratoolbox import behavioranalysis
import studyutils
import studyparams
from importlib import reload
reload(studyparams)
reload(studyutils)


DEBUG = 0
PLOT = 0
SAVE = 1

# -- Check if trialSubset should be changed. Options: 'running', 'notrunning' --
if len(sys.argv)>1:
    trialSubset = sys.argv[1]
    if trialSubset=='all':
        trialSubset = ''
else:
    trialSubset = ''
if trialSubset not in ['', 'all', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

if len(sys.argv)>2:
    evokedPeriod = sys.argv[2]
else:
    evokedPeriod = ''
if evokedPeriod not in ['', 'onset', 'sustained', 'offset']:
    raise ValueError("evokedPeriod must be '', 'onset', 'sustained', 'offset'")

N_PRE_ODDBALL = 8  # Analyze only the last N standard trials before oddball

if evokedPeriod == 'onset':
    timeRangeEvoked = [0.015, 0.065]
elif evokedPeriod == 'sustained':
    timeRangeEvoked = [0.065, 0.115]
elif evokedPeriod == 'offset':
    timeRangeEvoked = [0.115, 0.165]
else:
    timeRangeEvoked = [0.015, 0.115]
timeRange = {'Full':     [-0.3, 0.45],
             'Baseline': [-0.2, 0],
             'Evoked':   timeRangeEvoked }
timeRangeDuration = {k:np.diff(timeRange[k])[0] for k in timeRange.keys()}

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_coords.h5')
celldb = celldatabase.load_hdf(dbFilename)
#outputFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
if trialSubset == '':
    outputFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
else:
    outputFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_oddball_{trialSubset}.h5')
if evokedPeriod != '':
    outputFilename = outputFilename.replace('.h5', f'_{evokedPeriod}.h5')

# -- Load running data --
if trialSubset in ['running', 'notrunning']:
    runningData = studyutils.load_running_data()
    
# -- Initialize the dictionaries for new data (with keys like: doiDownStandardEvokedFiringRate) --
measurements = ['BaselineFiringRate', 'EvokedFiringRate', 'Ntrials']
stimConditions = ['oddball', 'standard']
columnsDict = {}
for reagent in studyparams.REAGENTS:
    for sessionType, sessionInfo in studyparams.ODDBALL_SESSION_TYPES.items():
        for stimCond in stimConditions:
            thisStim = sessionInfo[stimCond]
            columnsDict[reagent+thisStim+stimCond.capitalize()+'Pval'] = np.full(len(celldb), np.nan)
            for measurement in measurements:
                columnsDict[reagent+thisStim+stimCond.capitalize()+measurement] = np.full(len(celldb), np.nan)

if DEBUG:
    indRow = 46  # 46 # 55  # 1536
    celldbToUse = celldb.iloc[[indRow]]
else:
    celldbToUse = celldb
    
for indRow, dbRow in celldbToUse.iterrows():
    oneCell = ephyscore.Cell(dbRow)
    if PLOT:
        indplot = 0; plt.clf()
        plt.suptitle(f'{oneCell} [{indRow}]', fontweight='bold')
    if indRow%20==0:
        print(f'{indRow}/{len(celldb)} cells analyzed')
    for reagent in studyparams.REAGENTS:
        for sessionType, sessionInfo in studyparams.ODDBALL_SESSION_TYPES.items():
            ephysData, bdata = oneCell.load(reagent+sessionType)
            spikeTimes = ephysData['spikeTimes']
            eventOnsetTimes = ephysData['events']['stimOn']

            stimEachTrial = bdata['currentStartFreq']
            nTrials = len(stimEachTrial)

            # -- Doing this before selecting trials by running/notrunning --
            possibleStim = np.unique(stimEachTrial)
            nStim = len(possibleStim)

            # If the ephys data is 1 more than the bdata, delete the last ephys trial.
            if len(stimEachTrial) == len(eventOnsetTimes)-1:
                eventOnsetTimes = eventOnsetTimes[:len(stimEachTrial)]
            assert len(stimEachTrial) == len(eventOnsetTimes), \
                "Number of trials in behavior and ephys do not match for {oneCell}"
            
            # -- Select trials according to running state --
            behaviorDataFilePath = oneCell.get_behavior_path(reagent+sessionType)
            behaviorDataFilename = os.path.basename(behaviorDataFilePath)
            if trialSubset in ['running', 'notrunning']:
                if behaviorDataFilename not in runningData[oneCell.subject]:
                    continue
                runningThisSession = runningData[oneCell.subject][behaviorDataFilename]
                if trialSubset == 'running':
                    runningTrials = runningThisSession >= studyparams.RUNNING_THRESHOLD
                    eventOnsetTimes = eventOnsetTimes[runningTrials]
                    stimEachTrial = stimEachTrial[runningTrials]
                elif trialSubset == 'notrunning':
                    notRunningTrials = runningThisSession < studyparams.RUNNING_THRESHOLD
                    eventOnsetTimes = eventOnsetTimes[notRunningTrials]
                    stimEachTrial = stimEachTrial[notRunningTrials]
                nTrials = len(stimEachTrial)
                if nTrials==0:
                    continue
            
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange['Full'])

            trialsEachCond = behavioranalysis.find_trials_each_type(stimEachTrial, possibleStim)
            
            # -- Get oddball and standard trials (checking things match ODDBALL_SESSION_TYPES) --
            indOddball = np.argmin(trialsEachCond.sum(axis=0))
            indEachCond = {'oddball': indOddball, 'standard': 1-indOddball}
            # -- According to currentStartFreq, this is how the stim are ordered --
            stimOrder = {'FM': ['FM_up', 'FM_down'],
                         'Chord': ['low_freq', 'high_freq']}
            # -- Sanity check --
            oddballStim = bdata.labels['oddballStim'][bdata['oddballStim'][-1]]
            stimType = bdata.labels['stimType'][bdata['stimType'][-1]]
            assert oddballStim == stimOrder[stimType][indOddball], \
                f"Oddball stim should be {stimOrder[stimType][indOddball]}"

            # -- Include only the last N standard trials before oddball--
            oddballInds = np.where(trialsEachCond[:,indOddball])[0]
            zeroind, pre = np.meshgrid(oddballInds, np.arange(1, N_PRE_ODDBALL+1))
            lastXinds = np.concatenate(zeroind-pre)
            trialsEachCond[:, 1-indOddball] = False
            trialsEachCond[lastXinds, 1-indOddball] = True

            # -- Estimate firing rate --
            for stimcond in stimConditions:
                trialsThisCond = trialsEachCond[:, indEachCond[stimcond]]
                indexLimitsThisCond = indexLimitsEachTrial[:, trialsThisCond]
                firingRateDict = {}
                thisStim = sessionInfo[stimcond]
                thisKey = reagent + thisStim + stimcond.capitalize()
                nTrials = np.sum(trialsThisCond)
                columnsDict[thisKey+'Ntrials'][indRow] = nTrials
                for period in ['Baseline', 'Evoked']:
                    spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                             indexLimitsThisCond,
                                                                             timeRange[period])
                    firingRate = spikeCountMat.mean()/timeRangeDuration[period]
                    #thisStim = sessionInfo[stimcond]
                    #thisKey = reagent + thisStim + stimcond.capitalize() + period
                    columnsDict[thisKey+period+'FiringRate'][indRow] = firingRate
                    firingRateDict[period] = spikeCountMat.flatten()/timeRangeDuration[period]
                    #print(f'{thisKey}: \t{firingRate:0.3f} Hz')
                try:
                    wStat, pValThisCond = stats.wilcoxon(firingRateDict['Baseline'],
                                                         firingRateDict['Evoked'])
                except ValueError:
                    pValThisCond = np.nan
                #pValKey = reagent + thisStim + stimcond.capitalize() + 'Pval'
                columnsDict[thisKey + 'Pval'][indRow] = pValThisCond
                if DEBUG:
                    print(f'{thisKey}: p = \t{pValThisCond:0.3f}')
            
            if PLOT:
                try:
                    indplot += 1
                except:
                    #indplot = 1; plt.clf()
                    pass
                plt.subplot(3,4,indplot)
                pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                               timeRange['Full'], trialsEachCond)
                plt.setp(pRaster, ms=0.5)
                plt.title(f'{reagent} {sessionType}')
                plt.show()
                #plt.pause(0.5); 
    #plt.waitforbuttonpress()

celldbWithOddball = celldb.assign(**columnsDict)

# -- Save the updated celldb --
if SAVE:
    celldatabase.save_hdf(celldbWithOddball, outputFilename)

# -- Useful for debugging --    
# for k,v in celldbWithOddball.iloc[46].items(): print(f'{k}:\t {v}')
