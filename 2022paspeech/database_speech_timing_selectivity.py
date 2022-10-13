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
figDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'vot_timing')
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning_combineAudDAudPo.h5')
fulldb = celldatabase.load_hdf(dbPath)
newdbPath = os.path.join(databaseDir, 'fulldb_paspeech_latencies_ftvot.h5')
nCells = len(fulldb)

allSubjects = studyparams.EPHYS_MICE
#allSubjects = studyparams.TEST_MOUSE


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

        FTParamsEachTrial = bdata['targetFTpercent']
        possibleFTParams = np.unique(FTParamsEachTrial)
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        nFT = len(possibleFTParams)
        nVOT = len(possibleVOTParams)

        if len(VOTParamsEachTrial) == len(eventOnsetTimes)-1:
            eventOnsetTimes = eventOnsetTimes[:len(VOTParamsEachTrial)]
            print(f'[{indRow}] Warning! BehavTrials ({len(VOTParamsEachTrial)}) and ' +
                  f'EphysTrials ({len(eventOnsetTimes)})')

        #period = [0,0.12]
        #periodDuration = [x[1]-x[0] for x in allPeriods]
        #[periodName] = ['respOnset']
        #trialsEachCond = behavioranalysis.find_trials_each_type(VOTParamsEachTrial, possibleVOTParams)
        trialsEachCond = behavioranalysis.find_trials_each_combination(FTParamsEachTrial, possibleFTParams, VOTParamsEachTrial, possibleVOTParams)
        # calculate spike times
        #(spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \ spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)
        votLatencies_FTmin = np.ones(nVOT)
        votLatencies_FTmax = np.ones(nVOT)
        ftLatencies_VOTmin = np.ones(nVOT)
        ftLatencies_VOTmax = np.ones(nVOT)

        for indVOT, thisVOT in enumerate(possibleVOTParams):

            selected_VOTtrials_FTmin = trialsEachCond[:,0,indVOT]
            selected_VOTtrials_FTmax = trialsEachCond[:,3,indVOT]
            selected_FTtrials_VOTmin = trialsEachCond[:,indVOT,0]
            selected_FTtrials_VOTmax = trialsEachCond[:,indVOT,3]

            #VOT
            (spikeTimesFromEventOnsetVOT_FTmin, trialIndexForEachSpikeVOT_FTmin, indexLimitsEachTrialVOT_FTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmin], latencyTimeRange)
            (spikeTimesFromEventOnsetVOT_FTmax, trialIndexForEachSpikeVOT_FTmax, indexLimitsEachTrialVOT_FTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmax], latencyTimeRange)
            #FT
            (spikeTimesFromEventOnsetFT_VOTmin, trialIndexForEachSpikeFT_VOTmin, indexLimitsEachTrialFT_VOTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmin], latencyTimeRange)
            (spikeTimesFromEventOnsetFT_VOTmax, trialIndexForEachSpikeFT_VOTmax, indexLimitsEachTrialFT_VOTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmax], latencyTimeRange)

            # Estimate latency
            #VOT
            smoothWin = signal.windows.hann(7)
            try:
                respLatencyVOT_FTmin, interimVOT_FTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmin, indexLimitsEachTrialVOT_FTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                votLatencies_FTmin[indVOT] = respLatencyVOT_FTmin
            except:
                votLatencies_FTmin[indVOT] = -1
            try:
                respLatencyVOT_FTmax, interimVOT_FTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrialVOT_FTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                votLatencies_FTmax[indVOT] = respLatencyVOT_FTmax
            except:
                votLatencies_FTmax[indVOT] = -1

            #FT
            try:
                respLatencyFT_VOTmin, interimFT_VOTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmin, indexLimitsEachTrialFT_VOTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                ftLatencies_VOTmin[indVOT] = respLatencyFT_VOTmin
            except:
                ftLatencies_VOTmin[indVOT] = -1
            try:
                respLatencyFT_VOTmax, interimFT_VOTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmax, indexLimitsEachTrialFT_VOTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                ftLatencies_VOTmax[indVOT] = respLatencyFT_VOTmax
            except:
                votLatencies_FTmax[indVOT] = -1



            if indCell % 10 == 0:
                print(f'{indCell}/{len(celldbResp)}')
                pass
                #print(f'[{indRow}] {str(oneCell)}')



        ## FIGURES ##
        if 0: #votLatencies_FTmin[0]>0 and votLatencies_FTmin[-1]>0:
            plt.clf()
            gsMain = gs.GridSpec(2, 5)
            gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.1, wspace=0.4, hspace=0.4)
            markerSize = 2
            ax0 = plt.subplot(gsMain[0,0])
            selectedTrials_FTmin = trialsEachCond[:,0,:]

            (spikeTimesFromEventOnsetVOT_FTmin, trialIndexForEachSpikeVOT_FTmin, indexLimitsEachTrialVOT_FTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials_FTmin], latencyTimeRange)
            #plt.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=markerSize)
            nTrials = len(indexLimitsEachTrialVOT_FTmin)
            (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selectedTrials_FTmin, nTrials)
            lastTrialEachCond = np.cumsum(nTrialsEachCond)
            firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
            pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnsetVOT_FTmin, indexLimitsEachTrialVOT_FTmin, latencyTimeRange, selectedTrials_FTmin, labels = possibleVOTParams)
            plt.setp(pRaster, ms = 2)
            plt.axvline(0, color='b')
            plt.title('FTmin')
            #plt.axvline(respLatency, color='r')
            for indVOT, thisVOT in enumerate(possibleVOTParams):
                plt.plot([votLatencies_FTmin[indVOT], votLatencies_FTmin[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
                ax1 = plt.subplot(gsMain[0,5-indVOT], sharex=ax0)
                try:
                    respLatency, interim = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmin, indexLimitsEachTrialVOT_FTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                    plt.plot(interim['timeVec'], interim['psth'])
                except:
                    print('no spikes for this condition')

            ax0 = plt.subplot(gsMain[1,0])
            selectedTrials_FTmax = trialsEachCond[:,3,:]
            (spikeTimesFromEventOnsetVOT_FTmax, trialIndexForEachSpike_FTmax, indexLimitsEachTrial_FTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selectedTrials_FTmax], latencyTimeRange)
            nTrials = len(indexLimitsEachTrial_FTmax)
            (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selectedTrials_FTmax, nTrials)
            lastTrialEachCond = np.cumsum(nTrialsEachCond)
            firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
            pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrial_max, latencyTimeRange, selectedTrials_FTmax, labels = possibleVOTParams)
            plt.title('FTmax')
            plt.setp(pRaster, ms = 2)
            plt.axvline(0, color='b')
            #plt.axvline(respLatency, color='r')
            for indVOT, thisVOT in enumerate(possibleVOTParams):
                plt.plot([votLatencies_FTmax[indVOT], votLatencies_FTmax[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
                ax1 = plt.subplot(gsMain[1,5-indVOT], sharex=ax0)
                try:
                    respLatency, interim = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrial_FTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                    plt.plot(interim['timeVec'], interim['psth'])
                except:
                    print('no spikes for this condition')


            #plt.xlim([-0.2, 0.2])
            plt.suptitle(f'[{indRow}]: {1e3*respLatency:0.1f} ms')

            plt.show()
            #sys.exit()
            #plt.waitforbuttonpress()
            #sys.exit()
            #input("press enter for next cell")
            figName = f'{subject}_{dbRow.date}_{dbRow.maxDepth}um_c{dbRow.cluster}vot.png'
            figPath = os.path.join(figDir, figName)
            #plt.savefig(figPath, format='png')

            plt.close()
            #plt.pause(0.5)

        fulldb.at[indRow, 'respLatency_VOTmin_FTmin'] = votLatencies_FTmin[0]
        fulldb.at[indRow, 'respLatency_VOTmax_FTmin'] = votLatencies_FTmin[-1]
        fulldb.at[indRow, 'respLatency_VOTdiff_FTmin'] = (votLatencies_FTmin[-1] - votLatencies_FTmin[0])

        fulldb.at[indRow, 'respLatency_VOTmin_FTmax'] = votLatencies_FTmax[0]
        fulldb.at[indRow, 'respLatency_VOTmax_FTmax'] = votLatencies_FTmax[-1]
        fulldb.at[indRow, 'respLatency_VOTdiff_FTmax'] = (votLatencies_FTmax[-1] - votLatencies_FTmax[0])

        fulldb.at[indRow, 'respLatency_FTmin_VOTmin'] = ftLatencies_VOTmin[0]
        fulldb.at[indRow, 'respLatency_FTmax_VOTmin'] = ftLatencies_VOTmin[-1]
        fulldb.at[indRow, 'respLatency_FTdiff_VOTmin'] = (ftLatencies_VOTmin[-1] - ftLatencies_VOTmax[0])

        fulldb.at[indRow, 'respLatency_FTmin_VOTmax'] = ftLatencies_VOTmax[0]
        fulldb.at[indRow, 'respLatency_FTmax_VOTmax'] = ftLatencies_VOTmax[-1]
        fulldb.at[indRow, 'respLatency_FTdiff_VOTmax'] = (ftLatencies_VOTmax[-1] - ftLatencies_VOTmin[0])
celldatabase.save_hdf(fulldb, newdbPath)
