import os
import sys
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from jaratoolbox import colorpalette as cp
from scipy import signal
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from jaratoolbox import extraplots
sys.path.append('..')
import studyparams
import figparams

plt.rcParams['font.family'] = 'Helvetica'
plt.rcParams['svg.fonttype'] = 'none'  # To render as font rather than outlines
SAVE_FIGURE = 1
LATENCY_LINES = 0
#figDir = 'C:\\Users\\jenny\\tmp\\rasters'
figDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'timing_rasters')
figFormat = 'png' # 'pdf' or 'svg'
figSize = [14, 12]

## load database
databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning_allcells.h5')
fulldb = celldatabase.load_hdf(dbPath)
allSubjects = studyparams.EPHYS_MICE
#allSubjects = studyparams.TEST_MOUSE

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse
    celldb = fulldb[fulldb['subject']==thisMouse]
    nCells = len(celldb)

    inforecFile = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')

    for indRow, dbRow in celldb.iterrows():
        plt.clf()

        oneCell = ephyscore.Cell(dbRow)
        gsMain = gs.GridSpec(5, 4, height_ratios=(5, 1, 1, 1, 1))
        gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.1, wspace=0.4, hspace=0.5) #Change spacing of things
        plt.gcf().set_size_inches([14, 12])

        bestFtIndex = np.round(dbRow[['selectivityIndexFT_VOTmin', 'selectivityIndexFT_VOTmax']].values.max(0),2)
        bestVotIndex = np.round(dbRow[['selectivityIndexVOT_FTmin', 'selectivityIndexVOT_FTmin']].values.max(0),2)

        plt.suptitle(f'Recording Site:{dbRow.recordingSiteName}\n' + f'VOT-SI ={bestVotIndex}'+ f'FT-SI = {bestFtIndex},', fontsize=14, fontweight='bold', y = 0.99)
        VOTlabels = ['0', '20', '40', '60']
        FTlabels = ['9', '3', '-3', '-9']
        colorsEachVOT = [cp.TangoPalette['ScarletRed3'], cp.TangoPalette['ScarletRed2'], cp.TangoPalette['Butter2'], cp.TangoPalette['Butter3']]
        colorsEachFT = [cp.TangoPalette['SkyBlue3'], cp.TangoPalette['SkyBlue2'], cp.TangoPalette['Chameleon2'], cp.TangoPalette['Chameleon3']]
        #FTVOTBorders
        ephysData, bdata = oneCell.load('FTVOTBorders')

        # Align spikes to an event
        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']
        timeRange = [-0.3, 0.45]  # In seconds
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

        # Type-sorted rasters -- FTVOTBorders
        timeRange = [-0.075, 0.15]
        VOTParamsEachTrial = bdata['targetVOTpercent']
        possibleVOTParams = np.unique(VOTParamsEachTrial)
        FTParamsEachTrial = bdata['targetFTpercent']
        possibleFTParams = np.unique(FTParamsEachTrial)
        trialsEachCond = behavioranalysis.find_trials_each_combination(VOTParamsEachTrial, possibleVOTParams, FTParamsEachTrial, possibleFTParams)
        trialsEachVOT_FTmin = trialsEachCond[:, :, 0]
        trialsEachVOT_FTmax = trialsEachCond[:, :, -1]
        trialsEachFT_VOTmin = trialsEachCond[:, 0, :]
        trialsEachFT_VOTmax = trialsEachCond[:, -1, :]
        nVOT = len(possibleVOTParams)

        # Raster -- VOT (FTmin)
        ax0 = plt.subplot(gsMain[0, 0])
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmin, colorsEachVOT, labels = VOTlabels)
        plt.setp(pRaster, ms=2)
        #plt.xlabel('Time (s)')
        plt.ylabel('VOT (ms)')
        #plt.title(f'VOT, FTmin (n = {np.sum(selected_VOTtrials_FTmin)})')
        plt.title(r'$\Delta$VOT, FTmin')

        # Raster -- VOT (FTmax)
        ax1 = plt.subplot(gsMain[0, 1])
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmax, colorsEachVOT, labels = VOTlabels)
        plt.setp(pRaster, ms=2)
        #plt.xlabel('Time (s)')
        plt.ylabel('VOT (ms)')
        plt.title(r'$\Delta$VOT, FTmax')

        # Raster -- FT (VOT = min)
        ax2 = plt.subplot(gsMain[0, 2])
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmin, colorsEachFT, labels = FTlabels)
        plt.setp(pRaster, ms=2)
        #plt.xlabel('Time (s)')
        plt.ylabel('FT slope (oct/s)')
        plt.title(r'$\Delta$FT, VOTmin')


        # Raster -- FT (VOT = max)
        ax3 = plt.subplot(gsMain[0, 3])
        pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmax, colorsEachFT, labels = FTlabels)
        plt.setp(pRaster, ms=2)
        #plt.xlabel('Time (s)')
        plt.ylabel('FT slope (oct/s)')
        plt.title(r'$\Delta$FT, VOTmax')

        # PSTHs
        votLatencies_FTmin = np.ones(nVOT)
        votLatencies_FTmax = np.ones(nVOT)
        ftLatencies_VOTmin = np.ones(nVOT)
        ftLatencies_VOTmax = np.ones(nVOT)
        latencyTimeRange = [-0.15, 0.15]

        for indVOT, thisVOT in enumerate(possibleVOTParams):

            selected_VOTtrials_FTmin = trialsEachCond[:, indVOT, 0]
            selected_VOTtrials_FTmax = trialsEachCond[:, indVOT, 3]
            selected_FTtrials_VOTmin = trialsEachCond[:, 0, indVOT]
            selected_FTtrials_VOTmax = trialsEachCond[:, 3, indVOT]

            #VOT
            (spikeTimesFromEventOnsetVOT_FTmin, trialIndexForEachSpikeVOT_FTmin, indexLimitsEachTrialVOT_FTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmin], latencyTimeRange)
            (spikeTimesFromEventOnsetVOT_FTmax, trialIndexForEachSpikeVOT_FTmax, indexLimitsEachTrialVOT_FTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmax], latencyTimeRange)
            #FT
            (spikeTimesFromEventOnsetFT_VOTmin, trialIndexForEachSpikeFT_VOTmin, indexLimitsEachTrialFT_VOTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmin], latencyTimeRange)
            (spikeTimesFromEventOnsetFT_VOTmax, trialIndexForEachSpikeFT_VOTmax, indexLimitsEachTrialFT_VOTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmax], latencyTimeRange)
            plt.show()

            # ---Estimate latency---
            #--VOT--
            smoothWin = signal.windows.hann(13)
            yMax = np.round(dbRow[['speechFiringRateMaxOnset', 'speechFiringRateMaxSustain']].values.max(0))*8
            if yMax < 5:
                yMax = 10
            try:
                respLatencyVOT_FTmin, interimVOT_FTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmin, indexLimitsEachTrialVOT_FTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                votLatencies_FTmin[indVOT] = respLatencyVOT_FTmin
                #whichAx = f'ax{4*(indVOT+1)}'
                plt.subplot(gsMain[4-indVOT,0], sharex=ax0)
                timeBin = interimVOT_FTmin['timeVec'][1]-interimVOT_FTmin['timeVec'][0]
                psth = interimVOT_FTmin['psth']/timeBin
                plt.plot(interimVOT_FTmin['timeVec'], psth, color=colorsEachVOT[indVOT])
                plt.title(f'VOT {VOTlabels[indVOT]}ms')

                plt.ylim((0,yMax))
                if indVOT == 0:
                    plt.xlabel('Time (s)')
                if LATENCY_LINES:
                    nTrials = len(indexLimitsEachTrialVOT_FTmin[1,:])
                    (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_VOTtrials_FTmin, nTrials)
                    lastTrialEachCond = np.cumsum(nTrialsEachCond)
                    firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                    ax0.plot([votLatencies_FTmin[indVOT], votLatencies_FTmin[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
            except:
                votLatencies_FTmin[indVOT] = -1
            try:
                respLatencyVOT_FTmax, interimVOT_FTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrialVOT_FTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                votLatencies_FTmax[indVOT] = respLatencyVOT_FTmax
                plt.subplot(gsMain[4-indVOT,1], sharex=ax1)
                timeBin = interimVOT_FTmax['timeVec'][1]-interimVOT_FTmax['timeVec'][0]
                psth = interimVOT_FTmax['psth']/timeBin
                plt.plot(interimVOT_FTmax['timeVec'], psth, color=colorsEachVOT[indVOT])
                plt.title(f'VOT {VOTlabels[indVOT]}ms')
                plt.ylim((0,yMax))
                if indVOT == 0:
                    plt.xlabel('Time (s)')
                if LATENCY_LINES:
                    (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_VOTtrials_FTmax, nTrials)
                    lastTrialEachCond = np.cumsum(nTrialsEachCond)
                    firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                    ax1.plot([votLatencies_FTmax[indVOT], votLatencies_FTmax[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])

            except:
                votLatencies_FTmax[indVOT] = -1

            #--FT--
            try:
                respLatencyFT_VOTmin, interimFT_VOTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmin, indexLimitsEachTrialFT_VOTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                ftLatencies_VOTmin[indVOT] = respLatencyFT_VOTmin
                plt.subplot(gsMain[4-indVOT,2], sharex=ax2)
                timeBin = interimFT_VOTmin['timeVec'][1]-interimFT_VOTmin['timeVec'][0]
                psth = interimFT_VOTmin['psth']/timeBin
                plt.plot(interimFT_VOTmin['timeVec'], psth, color=colorsEachFT[indVOT])
                plt.title(f'FT {FTlabels[indVOT]}oct/s')
                plt.ylim((0,yMax))
                if indVOT == 0:
                    plt.xlabel('Time (s)')
                if LATENCY_LINES:
                    (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_FTtrials_VOTmin, nTrials)
                    lastTrialEachCond = np.cumsum(nTrialsEachCond)
                    firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                    ax2.plot([ftLatencies_VOTmin[indVOT], ftLatencies_VOTmin[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
            except:
                ftLatencies_VOTmin[indVOT] = -1
            try:
                respLatencyFT_VOTmax, interimFT_VOTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmax, indexLimitsEachTrialFT_VOTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                ftLatencies_VOTmax[indVOT] = respLatencyFT_VOTmax
                plt.subplot(gsMain[4-indVOT,3], sharex=ax3)
                timeBin = interimFT_VOTmax['timeVec'][1]-interimFT_VOTmax['timeVec'][0]
                psth = interimFT_VOTmax['psth']/timeBin
                plt.plot(interimFT_VOTmax['timeVec'], psth, color=colorsEachFT[indVOT])
                plt.title(f'FT {FTlabels[indVOT]} oct/s')
                plt.ylim((0,yMax))
                if indVOT == 0:
                    plt.xlabel('Time (s)')
                if LATENCY_LINES:
                    (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_FTtrials_VOTmax, nTrials)
                    lastTrialEachCond = np.cumsum(nTrialsEachCond)
                    firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                    ax3.plot([ftLatencies_VOTmax[indVOT], ftLatencies_VOTmax[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
            except:
                ftLatencies_VOTmax[indVOT] = -1



        #plt.show()
        #input("press enter for next cell")

        if SAVE_FIGURE:
            #figPath = os.path.join(figDir, f'{subject}_{dbRow.date}_{dbRow.maxDepth}um_c{dbRow.cluster}_rasters.png')
            FIGNAME = f'{subject}_{dbRow.date}_{dbRow.maxDepth}um_c{dbRow.cluster}_rasters'
            #plt.savefig(figPath, format='png')
            extraplots.save_figure(FIGNAME, figFormat, figSize, figDir)

        #plt.close()
