'''
This script assesses timing-based selectivity to VOT
'''
import numpy as np
import os
import sys
import studyparams
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
from scipy import signal
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from jaratoolbox import extraplots
import studyutils
from importlib import reload
reload(studyutils)

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning.h5')
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning_combineAudDAudPo.h5')
fulldb = celldatabase.load_hdf(dbPath)
newdbPath = os.path.join(databaseDir, 'fulldb_paspeech_latencies.h5')
nCells = len(fulldb)

allSubjects = studyparams.EPHYS_MICE
#allSubjects = studyparams.TEST_MOUSE

fulldb['respLatency_VOTmax'] = np.nan
fulldb['respLatency_VOTmin'] = np.nan
fulldb['respLatency_VOTdiff'] = np.nan

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse
    #dbPath = os.path.join(databaseDir, f'{subject}_paspeech_speech_pval.h5')
    #celldb = celldatabase.load_hdf(dbPath)
    celldb = fulldb[fulldb['subject']==thisMouse]
    nCells = len(celldb)

    N_FT = 4 # HARDCODED
    N_VOT = 4 #HARDCODED
    N_SPEECH = 12 #HARDCODED

    correctedAlpha = 0.05/N_SPEECH
    celldbResp = celldb[(celldb.speechMinPvalOnset < correctedAlpha) | (celldb.speechMinPvalSustain < correctedAlpha)]

    indCell = -1
    for indRow, dbRow in celldbResp.iterrows():
        indCell += 1
        oneCell = ephyscore.Cell(dbRow)

        ephysData, bdata = oneCell.load('FTVOTBorders')

        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        timeRange = [-0.4, 0.55]  # In seconds
        latencyTimeRange = [-0.12, 0.12]

        #FTParamsEachTrial = bdata['targetFTpercent']
        #possibleFTParams = np.unique(FTParamsEachTrial)
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        #nFT = len(possibleFTParams)
        nVOT = len(possibleVOTParams)

        if len(VOTParamsEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(VOTParamsEachTrial)]
            print(f'[{indRow}] Warning! BehavTrials ({len(VOTParamsEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')

        #period = [0,0.12]
        #periodDuration = [x[1]-x[0] for x in allPeriods]
        #[periodName] = ['respOnset']
        trialsEachCond = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)
        #trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)
        # calculate spike times
        #(spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \ spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)
        votLatencies = np.ones(nVOT)
        for indVOT, thisVOT in enumerate (possibleVOTParams):

            selectedTrials = trialsEachCond[:,indVOT]

            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials], latencyTimeRange)

            # Estimate latency
            smoothWin = signal.windows.hann(7)
            respLatency, interim = spikesanalysis.response_latency(spikeTimesFromEventOnset, indexLimitsEachTrial, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
            votLatencies[indVOT] = respLatency
            if indCell % 10 == 0:
                print(f'{indCell}/{len(celldbResp)}')
                pass
            #print(f'[{indRow}] {str(oneCell)}')

        if 0:
            plt.clf()
            markerSize = 2
            ax0 = plt.subplot(2,1,1)
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, latencyTimeRange)
            #plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=markerSize)
            nTrials = len(indexLimitsEachTrial)
            (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(trialsEachCond, nTrials)
            lastTrialEachCond = np.cumsum(nTrialsEachCond)
            firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
            pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, latencyTimeRange, trialsEachCond, labels = possibleVOTParams)
            plt.setp(pRaster, ms = 2)
            plt.axvline(0, color='b')
            #plt.axvline(respLatency, color='r')
            for indVOT, thisVOT in enumerate(possibleVOTParams):
                plt.plot([votLatencies[indVOT], votLatencies[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
            #plt.xlim([-0.2, 0.2])
            plt.title(f'[{indRow}]: {1e3*respLatency:0.1f} ms')
            ax1 = plt.subplot(2,1,2, sharex=ax0)
            plt.plot(interim['timeVec'], interim['psth'])
            plt.show()
            #sys.exit()
            #plt.waitforbuttonpress()
            #sys.exit()
            input("press enter for next cell")

            plt.close()
            #plt.pause(0.5)

        fulldb.at[indRow, 'respLatency_VOTmin'] = votLatencies[0]
        fulldb.at[indRow, 'respLatency_VOTmax'] = votLatencies[-1]
        fulldb.at[indRow, 'respLatency_VOTdiff'] = (votLatencies[-1] - votLatencies[0])
celldatabase.save_hdf(fulldb, newdbPath)
