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
SAVE_FIGURE = 0
LATENCY_LINES = 0
figDir = 'C:\\Users\\jenny\\tmp\\rasters'
#figDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'timing_rasters','exampleCells')
figFormat = 'png' # 'pdf' or 'svg'
figSize = [14, 12]
fontSizeTitles = figparams.fontSizeTitles
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks

## load database
databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning_allcells.h5')
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning.h5')
fulldb = celldatabase.load_hdf(dbPath)
FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
shuffledDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle=True)

pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
selectivityIndexVOT_FTmax = figData['selectivityIndexVOT_FTmax']
selectivityIndexVOT_FTmin = figData['selectivityIndexVOT_FTmin']
selectivityIndexFT_VOTmax = figData['selectivityIndexFT_VOTmax']
selectivityIndexFT_VOTmin = figData['selectivityIndexFT_VOTmin']

speechResponsive = figData['speechResponsive']
excludeCells = figData['excludeCells']

#allSubjects = studyparams.EPHYS_MICE
allSubjects = studyparams.TEST_MOUSE
#exampleCells = [1155, 1330, 777, 162]
indCell = -1

for indMouse, thisMouse in enumerate(allSubjects):
    subject = thisMouse
    celldb = fulldb[fulldb['subject']==thisMouse]
    nCells = len(celldb)


    for indRow, dbRow in celldb.iterrows():

        indCell += 1
        print(f'indCell = {indCell}')
        if speechResponsive[indCell] &~ excludeCells[indCell]:
            #dbRow = celldb.loc[thisCell]
            oneCell = ephyscore.Cell(dbRow)
            subject = thisMouse
            inforecFile = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')

            plt.clf()
            gsMain = gs.GridSpec(5, 4, height_ratios=(5, 1, 1, 1, 1))
            #gsMain = gs.GridSpec(2, 4, height_ratios=(8, 2))
            gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.1, wspace=0.4, hspace=0.5) #Change spacing of things
            plt.gcf().set_size_inches([14, 12])

            bestFtIndex = figData['bestSelectivityIndexFt'][indCell]
            bestVotIndex = figData['bestSelectivityIndexVot'][indCell]

            plt.suptitle(f'{indCell}_{dbRow.subject}_{dbRow.date}_c{dbRow.cluster}\n'+ f'Recording Site:{dbRow.recordingSiteName}', fontsize=fontSizeTitles, fontweight='bold', y = 0.99)
            VOTlabels = ['0', '20', '40', '60']
            FTlabels = ['9', '3', '-3', '-9']
            colorsEachVOT = [cp.TangoPalette['ScarletRed3'], cp.TangoPalette['ScarletRed2'], cp.TangoPalette['Butter2'], cp.TangoPalette['Butter3']]
            colorsEachFT = [cp.TangoPalette['SkyBlue3'], cp.TangoPalette['SkyBlue2'], cp.TangoPalette['Chameleon2'], cp.TangoPalette['Chameleon3']]


            #--Load data FTVOTBorders
            ephysData, bdata = oneCell.load('FTVOTBorders')

            # Align spikes to an event
            spikeTimes = ephysData['spikeTimes']
            eventOnsetTimes = ephysData['events']['stimOn']
            FTParamsEachTrial = bdata['targetFTpercent']


            if (len(FTParamsEachTrial)>len(eventOnsetTimes)) or (len(FTParamsEachTrial)<len(eventOnsetTimes)-1):
                print(f'[{indRow}] Warning! BevahTrials ({len(rateEachTrial)}) and ' +
                      f'EphysTrials ({len(eventOnsetTimes)})')
                continue
            if len(FTParamsEachTrial) == len(eventOnsetTimes)-1:
                eventOnsetTimes = eventOnsetTimes[:len(FTParamsEachTrial)]

            timeRange = [-0.3, 0.45]  # In seconds
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

            # Type-sorted rasters -- FTVOTBorders
            timeRange = [-0.075, 0.15]
            VOTParamsEachTrial = bdata['targetVOTpercent']
            possibleVOTParams = np.unique(VOTParamsEachTrial)
            possibleFTParams = np.unique(FTParamsEachTrial)


            trialsEachCond = behavioranalysis.find_trials_each_combination(VOTParamsEachTrial, possibleVOTParams, FTParamsEachTrial, possibleFTParams)
            trialsEachVOT_FTmin = trialsEachCond[:, :, 0]
            trialsEachVOT_FTmax = trialsEachCond[:, :, -1]
            trialsEachFT_VOTmin = trialsEachCond[:, 0, :]
            trialsEachFT_VOTmax = trialsEachCond[:, -1, :]
            nVOT = len(possibleVOTParams)
            nFT = len(possibleFTParams)
            pointSize = 4

            # Raster -- VOT (FTmin)
            ax0 = plt.subplot(gsMain[0, 0])
            pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmin, colorsEachVOT, labels = VOTlabels)
            plt.setp(pRaster, ms=pointSize)
            plt.xticks([])
            #plt.xlabel('Time (s)')
            plt.ylabel('VOT (ms)', fontsize=fontSizeLabels, fontweight='bold')
            #plt.title(f'VOT, FTmin (n = {np.sum(selected_VOTtrials_FTmin)})')
            plt.title(rf'$\Delta$VOT, FTmin. SI={np.round(selectivityIndexVOT_FTmin[indCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')

            # Raster -- VOT (FTmax)
            ax1 = plt.subplot(gsMain[0, 1])
            pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachVOT_FTmax, colorsEachVOT, labels = VOTlabels)
            plt.setp(pRaster, ms=pointSize)
            plt.xticks([])
            #plt.xlabel('Time (s)')
            plt.ylabel('VOT (ms)', fontsize=20, fontweight='bold')
            plt.title(rf'$\Delta$VOT, FTmax. SI={np.round(selectivityIndexVOT_FTmax[indCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')

            # Raster -- FT (VOT = min)
            ax2 = plt.subplot(gsMain[0, 2])
            pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmin, colorsEachFT, labels = FTlabels)
            plt.setp(pRaster, ms=pointSize)
            plt.xticks([])
            #plt.xlabel('Time (s)')
            plt.ylabel('FT slope (oct/s)', fontsize=20, fontweight='bold')
            plt.title(rf'$\Delta$FT, VOTmin. SI={np.round(selectivityIndexFT_VOTmin[indCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')


            # Raster -- FT (VOT = max)
            ax3 = plt.subplot(gsMain[0, 3])
            pRaster, hcond, zline =extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange, trialsEachFT_VOTmax, colorsEachFT, labels = FTlabels)
            plt.setp(pRaster, ms=pointSize)
            plt.xticks([])
            #plt.xlabel('Time (s)')
            plt.ylabel('FT slope (oct/s)', fontsize=20, fontweight='bold')
            plt.title(rf'$\Delta$FT, VOTmax. SI={np.round(selectivityIndexFT_VOTmax[indCell], 2)}', fontsize=fontSizeTitles, fontweight='bold')

            # PSTHs
            votLatencies_FTmin = np.ones(nVOT)
            votLatencies_FTmax = np.ones(nVOT)
            ftLatencies_VOTmin = np.ones(nVOT)
            ftLatencies_VOTmax = np.ones(nVOT)
            latencyTimeRange = [-0.15, 0.15]
            period = [0, 0.12]
            periodDuration = [period[1] - period[0]]
            spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset, indexLimitsEachTrial, period)
            spikesEachTrial = spikeCountMat[:,0]
            meanFiringRateOnset = np.empty([nFT, nVOT])
            for indFT, thisFT in enumerate(possibleFTParams):
                for indVOT, thisVOT in enumerate(possibleVOTParams):
                    if ((indFT == 1)| (indFT==2)) & ((indVOT==1)|(indVOT==2)):
                        meanFiringRateOnset[indFT, indVOT] = -1
                    else:
                        trialsThisCond = trialsEachCond[:, indVOT, indFT]
                        firingRateOnset = spikesEachTrial[trialsThisCond]/periodDuration
                        meanFiringRateOnset[indFT, indVOT] = firingRateOnset.mean()

            for indVOT, thisVOT in enumerate(possibleVOTParams):

                selected_VOTtrials_FTmin = trialsEachCond[:, indVOT, 0]
                selected_VOTtrials_FTmax = trialsEachCond[:, indVOT, -1]
                selected_FTtrials_VOTmin = trialsEachCond[:, 0, indVOT]
                selected_FTtrials_VOTmax = trialsEachCond[:, -1, indVOT]

                #VOT
                (spikeTimesFromEventOnsetVOT_FTmin, trialIndexForEachSpikeVOT_FTmin, indexLimitsEachTrialVOT_FTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmin], latencyTimeRange)
                (spikeTimesFromEventOnsetVOT_FTmax, trialIndexForEachSpikeVOT_FTmax, indexLimitsEachTrialVOT_FTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_VOTtrials_FTmax], latencyTimeRange)
                #FT
                (spikeTimesFromEventOnsetFT_VOTmin, trialIndexForEachSpikeFT_VOTmin, indexLimitsEachTrialFT_VOTmin) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmin], latencyTimeRange)
                (spikeTimesFromEventOnsetFT_VOTmax, trialIndexForEachSpikeFT_VOTmax, indexLimitsEachTrialFT_VOTmax) = spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes[selected_FTtrials_VOTmax], latencyTimeRange)
                plt.show()

                # ---Estimate latency---
                #--VOT--
                smoothWin = signal.windows.hann(9)
                #yMax = np.round(dbRow[['speechFiringRateMaxOnset', 'speechFiringRateMaxSustain', 'speechFiringRateBaseline']].values.max(0))*14
                #if yMax < 5:
                #    yMax = 10
                #try:
                respLatencyVOT_FTmin, interimVOT_FTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmin, indexLimitsEachTrialVOT_FTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                respLatencyVOT_FTmax, interimVOT_FTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrialVOT_FTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                respLatencyFT_VOTmin, interimFT_VOTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmin, indexLimitsEachTrialFT_VOTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                respLatencyFT_VOTmax, interimFT_VOTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmax, indexLimitsEachTrialFT_VOTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)

                timeBin = interimVOT_FTmin['timeVec'][1]-interimVOT_FTmin['timeVec'][0]

                yMax = (np.max([interimVOT_FTmin['maxResponse'], interimVOT_FTmax['maxResponse'], interimFT_VOTmin['maxResponse'], interimFT_VOTmax['maxResponse']])/timeBin)*1.8

                votLatencies_FTmin[indVOT] = respLatencyVOT_FTmin
                votLatencies_FTmax[indVOT] = respLatencyVOT_FTmax
                ftLatencies_VOTmin[indVOT] = respLatencyFT_VOTmin
                ftLatencies_VOTmax[indVOT] = respLatencyFT_VOTmax



                #whichAx = f'ax{4*(indVOT+1)}'
                plt.subplot(gsMain[4-indVOT,0], sharex=ax0)
                #plt.subplot(gsMain[1,0], sharex=ax0)
                #timeBin = interimVOT_FTmin['timeVec'][1]-interimVOT_FTmin['timeVec'][0]
                psth = interimVOT_FTmin['psth']/timeBin
                plt.plot(interimVOT_FTmin['timeVec'], psth, color=colorsEachVOT[indVOT])
                plt.title(f'VOT {VOTlabels[indVOT]}ms', fontsize = fontSizeLabels)
                plt.annotate(f'FR = {np.round(meanFiringRateOnset[0, indVOT], 1)}', xy = (0.7, 0.8), xycoords = 'axes fraction')

                #plt.xlim(-0.1,.15)
                plt.ylim((0,yMax))
                #plt.xticks([])
                if indVOT == 0:
                    plt.xlabel('Time (ms)', fontsize = fontSizeLabels)
                if LATENCY_LINES:
                    nTrials = len(indexLimitsEachTrialVOT_FTmin[1,:])
                    (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_VOTtrials_FTmin, nTrials)
                    lastTrialEachCond = np.cumsum(nTrialsEachCond)
                    firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                    ax0.plot([votLatencies_FTmin[indVOT], votLatencies_FTmin[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
                #except:
                #    votLatencies_FTmin[indVOT] = -1
                try:
                    #respLatencyVOT_FTmax, interimVOT_FTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetVOT_FTmax, indexLimitsEachTrialVOT_FTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                    #votLatencies_FTmax[indVOT] = respLatencyVOT_FTmax
                    plt.subplot(gsMain[4-indVOT,1], sharex=ax1)
                    #timeBin = interimVOT_FTmax['timeVec'][1]-interimVOT_FTmax['timeVec'][0]
                    psth = interimVOT_FTmax['psth']/timeBin
                    plt.plot(interimVOT_FTmax['timeVec'], psth, color=colorsEachVOT[indVOT])
                    plt.title(f'VOT {VOTlabels[indVOT]}ms', fontsize = fontSizeLabels)
                    plt.annotate(f'FR = {np.round(meanFiringRateOnset[-1, indVOT], 1)}', xy = (0.7, 0.8), xycoords = 'axes fraction')
                    #plt.xlim(-0.1,.15)
                    plt.ylim((0,yMax))
                    #plt.xticks([])
                    if indVOT == 0:
                        plt.xlabel('Time (ms)', fontsize = fontSizeLabels)
                    if LATENCY_LINES:
                        (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_VOTtrials_FTmax, nTrials)
                        lastTrialEachCond = np.cumsum(nTrialsEachCond)
                        firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                        ax1.plot([votLatencies_FTmax[indVOT], votLatencies_FTmax[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])

                except:
                    votLatencies_FTmax[indVOT] = -1

                #--FT--
                try:
                    #respLatencyFT_VOTmin, interimFT_VOTmin = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmin, indexLimitsEachTrialFT_VOTmin, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                    #ftLatencies_VOTmin[indVOT] = respLatencyFT_VOTmin
                    plt.subplot(gsMain[4-indVOT,2], sharex=ax2)
                    #timeBin = interimFT_VOTmin['timeVec'][1]-interimFT_VOTmin['timeVec'][0]
                    psth = interimFT_VOTmin['psth']/timeBin
                    plt.plot(interimFT_VOTmin['timeVec'], psth, color=colorsEachFT[indVOT])
                    plt.title(f'FT {FTlabels[indVOT]}oct/s', fontsize = fontSizeTitles)
                    plt.annotate(f'FR = {np.round(meanFiringRateOnset[indVOT, 0],1)}', xy = (0.7, 0.8), xycoords = 'axes fraction')
                    #plt.xlim(-0.1,.15)
                    plt.ylim((0,yMax))
                    #plt.xticks([])
                    if indVOT == 0:
                        plt.xlabel('Time (ms)', fontsize = fontSizeLabels)
                    if LATENCY_LINES:
                        (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_FTtrials_VOTmin, nTrials)
                        lastTrialEachCond = np.cumsum(nTrialsEachCond)
                        firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                        ax2.plot([ftLatencies_VOTmin[indVOT], ftLatencies_VOTmin[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
                except:
                    ftLatencies_VOTmin[indVOT] = -1
                try:
                    #respLatencyFT_VOTmax, interimFT_VOTmax = spikesanalysis.response_latency(spikeTimesFromEventOnsetFT_VOTmax, indexLimitsEachTrialFT_VOTmax, latencyTimeRange, threshold = 0.5, win = smoothWin, invert = False)
                    #ftLatencies_VOTmax[indVOT] = respLatencyFT_VOTmax
                    plt.subplot(gsMain[4-indVOT,3], sharex=ax3)
                    #timeBin = interimFT_VOTmax['timeVec'][1]-interimFT_VOTmax['timeVec'][0]
                    psth = interimFT_VOTmax['psth']/timeBin
                    plt.plot(interimFT_VOTmax['timeVec'], psth, color=colorsEachFT[indVOT])
                    plt.title(f'FT {FTlabels[indVOT]} oct/s', fontsize = fontSizeLabels)
                    plt.annotate(f'FR = {np.round(meanFiringRateOnset[indVOT, -1],1)}', xy = (0.7, 0.8), xycoords = 'axes fraction')
                    plt.ylim((0,yMax))
                    #plt.xticks([])
                    if indVOT == 0:
                        plt.xlabel('Time (ms)', fontsize = fontSizeLabels)
                    if LATENCY_LINES:
                        (trialsEachCondInds, nTrialsEachCond, nCond) = extraplots.trials_each_cond_inds(selected_FTtrials_VOTmax, nTrials)
                        lastTrialEachCond = np.cumsum(nTrialsEachCond)
                        firstTrialEachCond = np.r_[0, lastTrialEachCond[:-1]]
                        ax3.plot([ftLatencies_VOTmax[indVOT], ftLatencies_VOTmax[indVOT]], [firstTrialEachCond[indVOT], lastTrialEachCond[indVOT]])
                except:
                    ftLatencies_VOTmax[indVOT] = -1

            #print(votLatencies_FTmax)
            #print(votLatencies_FTmin)

            plt.show()
            #input("press enter for next cell")

            if SAVE_FIGURE:
                #figPath = os.path.join(figDir, f'{subject}_{dbRow.date}_{dbRow.maxDepth}um_c{dbRow.cluster}_rasters.png')
                FIGNAME = f'{indCell}_{subject}_{dbRow.date}_{dbRow.maxDepth}um_c{dbRow.cluster}_rasters'
                #plt.savefig(figPath, format='png')
                extraplots.save_figure(FIGNAME, figFormat, figSize, figDir)

            #plt.close()
