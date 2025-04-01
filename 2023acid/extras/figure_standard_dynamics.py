"""
Plot the PSTH of each standard in the sequence (ending with the oddball).
"""

import sys
sys.path.append('..')
import os
import numpy as np
import matplotlib.pyplot as plt
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
from jaratoolbox import extraplots
from jaratoolbox import behavioranalysis
import studyparams
from importlib import reload
reload(studyparams)

timeRange = [-0.1, 0.3]
binWidth = 0.010
timeVec = np.arange(timeRange[0], timeRange[-1], binWidth)
smoothWinSizePsth = 2 
downsampleFactorPsth = 1
lwPsth = 2

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_coords.h5')
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Select cells --
firingRateThreshold = 20
selcells = ((celldb.preDownStandardEvokedFiringRate>firingRateThreshold) |
            (celldb.preUpStandardEvokedFiringRate>firingRateThreshold))
celldb = celldb[selcells]

stimConditions = ['oddball', 'standard']
reagentsToPlot = ['pre'] #studyparams.REAGENTS
keysToPlot = ['FM_Down', 'FM_Up']
#keysToPlot = ['HighFreq', 'LowFreq']
sessionTypesToPlot = {key: studyparams.ODDBALL_SESSION_TYPES[key] for key in keysToPlot}

def find_oddball_trials(trialsEachCond):
    """Return a list of indices of the oddball trials"""
    indOddball = np.argmin(trialsEachCond.sum(axis=0))
    oddballInds = np.flatnonzero(trialsEachCond[:, indOddball])
    return oddballInds


# 'acid006', '2023-05-18', 3000, g0c240  # Strong response.
print('Press the SPACEBAR for next cell. Cell 46 has a strong response')

#indRow = 46
#for dbRow in [celldb.iloc[indRow]]:
for indRow, dbRow in celldb.iterrows():
    plt.clf()
    oneCell = ephyscore.Cell(dbRow)
    for reagent in reagentsToPlot:
        for sessionType, sessionInfo in sessionTypesToPlot.items():
            ephysData, bdata = oneCell.load(reagent+sessionType)
            spikeTimes = ephysData['spikeTimes']
            eventOnsetTimes = ephysData['events']['stimOn']

            stimEachTrial = bdata['currentStartFreq']
            nTrials = len(stimEachTrial)
    
            # If the ephys data is 1 more than the bdata, delete the last ephys trial.
            if len(stimEachTrial) == len(eventOnsetTimes)-1:
                eventOnsetTimes = eventOnsetTimes[:len(stimEachTrial)]
            assert len(stimEachTrial) == len(eventOnsetTimes), \
                "Number of trials in behavior and ephys do not match for {oneCell}"
            
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

            possibleStim = np.unique(stimEachTrial)
            trialsEachCond = behavioranalysis.find_trials_each_type(stimEachTrial, possibleStim)

            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                     indexLimitsEachTrial, timeVec)

            
            # -- Get oddball and standard trials (checking things match ODDBALL_SESSION_TYPES) --
            #indOddball = np.argmin(trialsEachCond.sum(axis=0))
            #indEachCond = {'oddball': indOddball, 'standard': 1-indOddball}

            oddballInds = find_oddball_trials(trialsEachCond)

            for stimPosFromOdd in range(-8, 1):
                trialsThisStim = oddballInds+stimPosFromOdd
                if sessionType=='FM_Down' or sessionType=='HighFreq':
                    plt.subplot(2, 9, stimPosFromOdd+9)
                elif sessionType=='FM_Up' or sessionType=='LowFreq':
                    plt.subplot(2, 9, stimPosFromOdd+9+9)
                    #plt.setp(pPSTH, color='r')
                    #plt.setp(pPSTH, color='b')
                pPSTH = extraplots.plot_psth(spikeCountMat[trialsThisStim]/binWidth,
                                             smoothWinSizePsth, timeVec,
                                             linewidth=lwPsth, downsamplefactor=downsampleFactorPsth)
                plt.grid(True)
                plt.ylim([0, 130])
                #plt.title(f'{oneCell}')
    plt.suptitle(f'{oneCell} [{indRow}]')
    plt.draw()
    plt.waitforbuttonpress()

# bdata.events['eventTime'][bdata.events['nextState']==bdata.stateMatrix['statesNames']['output1On']]

# -- Time between one stim and the next is between 627ms and 705ms --
# 1000*np.max(np.diff(bdata.events['eventTime'][bdata.events['nextState']==bdata.stateMatrix['statesNames']['output1On']]))
# 1000*np.min(np.diff(bdata.events['eventTime'][bdata.events['nextState']==bdata.stateMatrix['statesNames']['output1On']]))
