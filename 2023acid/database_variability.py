"""
Estimate trial-to-trial variability.

NOTE: because I am loading celldb_..._freqtuning.h5, and calculating running/notrunning
      here, I'm not saving the correct number of trials in 'preToneNtrials'.
      Maybe I should load celldb_..._freqtuning_running.h5 instead.

When estimating values for running and notRunning conditions, some sessions have very
few trials, so we select only sessions that have at least 4 trials of each type.

* DEBUGGING:
- The following don't match because not all freq will have trials:
np.unique(bdata['currentFreq'][runningTrials], return_counts=True)
trialsEachCond.sum(axis=0)
- These ones match:
nSpikesEachStim
print(nSpikesThisStim)

This script adds columns to the database with variability results for each reagent.

To run for only running or not-running trials, use:
    python database_freq_tuning.py running
    python database_freq_tuning.py notrunning

In Santiago's desktop it took between <2.5 min to run.
"""

import os
import sys
import numpy as np
from scipy import stats
from scipy import optimize
import matplotlib.pyplot as plt
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
from jaratoolbox import extraplots
from jaratoolbox import behavioranalysis
import studyparams
import studyutils
from importlib import reload
reload(studyparams)
reload(spikesanalysis)
reload(ephyscore)


DEBUG = 0
PLOT = 0
SAVE = 1

# -- Check if trialSubset should be changed. Options: 'running', 'notrunning' --
if len(sys.argv)==2:
    trialSubset = sys.argv[1]
else:
    trialSubset = ''
if trialSubset not in ['', 'running', 'notrunning']:
    raise ValueError("trialSubset must be '', 'all', 'running', or 'notrunning'")

timeRange = {'Full':     [-0.3, 0.45],
             'Baseline': [-0.2, 0],
             'Evoked':   [0.015, 0.115] }
timeRangeDuration = {k:np.diff(timeRange[k])[0] for k in timeRange.keys()}

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldb = celldatabase.load_hdf(dbFilename)
if trialSubset == '':
    outputFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_variability.h5')
else:
    outputFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_variability_{trialSubset}.h5')

# -- Load running data --
if trialSubset in ['running', 'notrunning']:
    runningData = studyutils.load_running_data()
    
# -- Initialize the dictionaries for new data (with keys like: doiBaselineFiringRate) --
measurements = ['ToneTrialToTrialVarSpks', 'ToneTrialToTrialMeanSpks', 'TonePvalEvokedEachStim']
'''
measurements = ['ToneBaselineFiringRate', 'ToneResponseMinPval', 'ToneSelectivityPval',
                'ToneFiringRateBestFreq', 'ToneBestFreq', 'ToneAvgEvokedFiringRate',
                'ToneGaussianA', 'ToneGaussianX0', 'ToneGaussianSigma', 'ToneGaussianY0',
                'ToneGaussianRsquare', 'ToneNtrials']
'''
columnsDict = {}
for reagent in studyparams.REAGENTS:
    for measurement in measurements:
        columnsDict[reagent+measurement] = np.full((len(celldb),studyparams.N_FREQ), np.nan)
        #columnsDict[reagent+measurement] = np.full(len(celldb), np.nan)
    #columnsDict[reagent+'ToneFiringRateEachFreq'] = np.full((len(celldb),studyparams.N_FREQ), np.nan)

if DEBUG:
    indRow = 46  # 46 # 55
    #indRow = 53  # Inverted
    #indRow = 176
    #indRow = 1533
    #indRow = 1318 # U-shape
    #indRow = 1501 # Wide going down
    #indRow = 1583 # very flat
    #indRow = 67  # Did not fit before fixing p0(A)
    celldbToUse = celldb.iloc[[indRow]]
else:
    celldbToUse = celldb
    
for indRow, dbRow in celldbToUse.iterrows():
    oneCell = ephyscore.Cell(dbRow)
    if PLOT:
        indplot = 0;
        fig = plt.gcf()
        fig.clf()
        plt.suptitle(f'{oneCell} [{indRow}]', fontweight='bold')
        ax = fig.subplots(nrows=3, ncols=2)
    if indRow%20==0:
        print(f'{indRow}/{len(celldb)} cells analyzed')
    for indr, reagent in enumerate(studyparams.REAGENTS):
        sessionType = 'PureTones'

        ephysData, bdata = oneCell.load(reagent+sessionType)
        spikeTimes = ephysData['spikeTimes']
        eventOnsetTimes = ephysData['events']['stimOn']

        stimEachTrial = bdata['currentFreq']
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
                print(f'---------- {behaviorDataFilename} not in runningData ------------')
                continue
            runningThisSession = runningData[oneCell.subject][behaviorDataFilename]
            runningTrials = runningThisSession >= studyparams.RUNNING_THRESHOLD
            notRunningTrials = runningThisSession < studyparams.RUNNING_THRESHOLD
            # The code below should work even if we running are not the same as (not notRunning)
            countfun = behavioranalysis.find_trials_each_type
            nTrialsEachStimRunning = countfun(stimEachTrial[runningTrials],
                                              possibleStim).sum(axis=0)
            nTrialsEachStimNotRunning = countfun(stimEachTrial[notRunningTrials],
                                                 possibleStim).sum(axis=0)
            if trialSubset == 'running':
                eventOnsetTimes = eventOnsetTimes[runningTrials]
                stimEachTrial = stimEachTrial[runningTrials]
            elif trialSubset == 'notrunning':
                eventOnsetTimes = eventOnsetTimes[notRunningTrials]
                stimEachTrial = stimEachTrial[notRunningTrials]
            nTrials = len(stimEachTrial)
            if nTrials==0:
                continue
            
        (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
            spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange['Full'])

        trialsEachCond = behavioranalysis.find_trials_each_type(stimEachTrial, possibleStim)

        # -- Estimate baseline firing rate --
        spikeCountMatBase = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 timeRange['Baseline'])
        baselineFiringRate = spikeCountMatBase.mean()/timeRangeDuration['Baseline']

        # -- Estimate evoked firing rate for each stim --
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,
                                                                 indexLimitsEachTrial,
                                                                 timeRange['Evoked'])
        nSpikesEachTrial = spikeCountMat[:,0]  # Flatten it

        # -- Calculate nSpikes for each freq to test selectivity --
        nSpikesEachStim = []
        avgSpikesEachStim = np.empty(nStim)
        stdSpikesEachStim = np.empty(nStim)
        pValEvokedEachStim = np.empty(nStim)
        for indStim, frequency in enumerate(possibleStim):
            nSpikesThisStim = nSpikesEachTrial[trialsEachCond[:,indStim]]
            nSpikesEachStim.append(nSpikesThisStim)
            #avgFiringRateEachStim[indStim] = nSpikesThisStim.mean()/timeRangeDuration['Evoked']
            avgSpikesEachStim[indStim] = nSpikesThisStim.mean()
            stdSpikesEachStim[indStim] = nSpikesThisStim.std()
            # -- Calculate p-value for each stim --
            baselineSpikesThisStim = spikeCountMatBase[trialsEachCond[:,indStim],0]
            try:
                wStat, pValThisStim = stats.wilcoxon(nSpikesThisStim, baselineSpikesThisStim)
            except ValueError:
                pValThisStim = 1
            pValEvokedEachStim[indStim] = pValThisStim
        '''
        try:
            nTrialsEachCond = trialsEachCond.sum(axis=0)>0
            nSpks = [nSpikesEachStim[ind] for ind in np.flatnonzero(nTrialsEachCond)]
            kStat, pValKruskal = stats.kruskal(*nSpks)
            #kStat, pValKruskal = stats.kruskal(*nSpikesEachStim)
        except ValueError:
            pValKruskal = 1
        '''
        pValEvokedMin = np.min(pValEvokedEachStim)

        freqEachTrial = stimEachTrial
        possibleFreq = possibleStim
        logFreq = np.log2(freqEachTrial)
        possibleLogFreq = np.log2(possibleFreq)
        '''
        # -- Fit Gaussian to tuning data --
        freqEachTrial = stimEachTrial
        possibleFreq = possibleStim
        logFreq = np.log2(freqEachTrial)
        possibleLogFreq = np.log2(possibleFreq)
        maxEvokedFiringRate = np.nanmax(avgFiringRateEachStim)
        changeFromBaseline = avgFiringRateEachStim - baselineFiringRate
        maxChangeFromBaseline = changeFromBaseline[np.nanargmax(np.abs(changeFromBaseline))]
        avgEvokedFiringRate = np.nanmean(avgFiringRateEachStim)
        
        #baselineFiringRate
        # PARAMS: a, x0, sigma, y0
        minS = studyparams.MIN_SIGMA
        maxS = studyparams.MAX_SIGMA
        if maxChangeFromBaseline >= 0:
            p0 = [maxChangeFromBaseline, possibleLogFreq[nStim//2], 1, baselineFiringRate]
            bounds = ([0, possibleLogFreq[0], minS, 0],
                      [np.inf, possibleLogFreq[-1], maxS, np.inf])
        else:
            p0 = [maxChangeFromBaseline, possibleLogFreq[nStim//2], 1, baselineFiringRate]
            bounds = ([-np.inf, possibleLogFreq[0], minS, 0],
                      [0, possibleLogFreq[-1], maxS, np.inf])
        try:
            firingRateEachTrial = nSpikesEachTrial/timeRangeDuration['Evoked']
            fitParams, pcov = optimize.curve_fit(spikesanalysis.gaussian, logFreq,
                                                 firingRateEachTrial, p0=p0, bounds=bounds)
        except RuntimeError:
            print("Could not fit gaussian curve to tuning data.")
            fitParams = np.full(4, np.nan)
            Rsquared = 0
        else:
            gaussianResp = spikesanalysis.gaussian(logFreq, *fitParams)
            residuals = firingRateEachTrial - gaussianResp
            ssquared = np.sum(residuals**2)
            ssTotal = np.sum((firingRateEachTrial-np.mean(firingRateEachTrial))**2)
            Rsquared = 1 - (ssquared/ssTotal)
            fullWidthHalfMax = 2.355*fitParams[2] # Sigma is fitParams[2]

        # -- Store results in dictionary --
        columnsDict[reagent+'ToneBaselineFiringRate'][indRow] = baselineFiringRate
        columnsDict[reagent+'ToneResponseMinPval'][indRow] = pValEvokedMin
        columnsDict[reagent+'ToneSelectivityPval'][indRow] = pValKruskal
        columnsDict[reagent+'ToneFiringRateBestFreq'][indRow] = maxEvokedFiringRate
        columnsDict[reagent+'ToneAvgEvokedFiringRate'][indRow] = avgEvokedFiringRate
        columnsDict[reagent+'ToneBestFreq'][indRow] = possibleFreq[avgFiringRateEachStim.argmax()]
        columnsDict[reagent+'ToneGaussianA'][indRow] = fitParams[0]
        columnsDict[reagent+'ToneGaussianX0'][indRow] = fitParams[1]
        columnsDict[reagent+'ToneGaussianSigma'][indRow] = fitParams[2]
        columnsDict[reagent+'ToneGaussianY0'][indRow] = fitParams[3]
        columnsDict[reagent+'ToneGaussianRsquare'][indRow] = Rsquared
        columnsDict[reagent+'ToneNtrials'][indRow] = nTrials
        columnsDict[reagent+'ToneFiringRateEachFreq'][indRow] = avgFiringRateEachStim
        '''

        # Exclude stimuli with too few trials
        #_, nTrialsEachStimRunning = np.unique(bdata['currentFreq'][runningTrials], return_counts=True)
        #_, nTrialsEachStimNotRunning = np.unique(bdata['currentFreq'][notRunningTrials], return_counts=True)
        if trialSubset in ['running', 'notrunning']:
            minNtrials = 4
            tooFewTrials = ((nTrialsEachStimRunning<minNtrials) |
                            (nTrialsEachStimNotRunning<minNtrials))
            avgSpikesEachStim[tooFewTrials] = np.nan
            stdSpikesEachStim[tooFewTrials] = np.nan
        columnsDict[reagent+'ToneTrialToTrialMeanSpks'][indRow] = avgSpikesEachStim
        columnsDict[reagent+'ToneTrialToTrialVarSpks'][indRow] = stdSpikesEachStim**2
        columnsDict[reagent+'TonePvalEvokedEachStim'][indRow] = pValEvokedEachStim
        fanoFactor = (stdSpikesEachStim**2)/avgSpikesEachStim
        
        if PLOT: # and pValEvokedMin<0.001:
            #plt.suptitle(f'[{indRow}]  {oneCell}')
            #plt.subplot(1,2,1)
            plt.sca(ax[indr, 0])
            pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset,
                                                           indexLimitsEachTrial,
                                                           timeRange['Full'], trialsEachCond)
            plt.setp(pRaster, ms=0.5)
            plt.title(f'{reagent} {sessionType}')
            
            #plt.subplot(1,2,2)
            plt.sca(ax[indr, 1])            
            pdots = plt.plot(possibleLogFreq, avgSpikesEachStim, 'o')
            '''
            if not np.isnan(fitParams[0]):
                xvals = np.linspace(possibleLogFreq[0], possibleLogFreq[-1], 60)
                yvals = spikesanalysis.gaussian(xvals, *fitParams)
                pfit = plt.plot(xvals, yvals, '-', lw=3)
            '''
            plt.ylabel('Avg spike count')
            plt.xlabel('Frequency (kHz)')
            xTickLabels = [f'{freq/1000:0.0f}' for freq in possibleFreq]
            plt.xticks(possibleLogFreq, xTickLabels)
            #plt.title(f'R2={Rsquared:0.4f}  s={fitParams[2]:0.4f}')

            axVar = ax[indr, 1].twinx()
            
            axVar.plot(possibleLogFreq, fanoFactor, 'ro')
            axVar.set_ylabel('Fano factor', color='r')
            axVar.set_ylim([0, 1.5])
            avgFanoFactor = np.nanmean(fanoFactor)
            axVar.text(0.95, 0.9, f'FF: {avgFanoFactor:0.2f}', color='r',
                       transform=axVar.transAxes, ha='right', fontsize=10)
            plt.show()
            plt.pause(0.5);
            #print(fitParams)
            if reagent=='doi':
                sys.exit()

for reagent in studyparams.REAGENTS:
    columnsDict[reagent+'ToneTrialToTrialMeanSpks'] = list(columnsDict[reagent+'ToneTrialToTrialMeanSpks'])
    columnsDict[reagent+'ToneTrialToTrialVarSpks'] = list(columnsDict[reagent+'ToneTrialToTrialVarSpks'])
    columnsDict[reagent+'TonePvalEvokedEachStim'] = list(columnsDict[reagent+'TonePvalEvokedEachStim'])
celldbWithVariability = celldb.assign(**columnsDict)

# -- Save the updated celldb --
if SAVE:
    celldatabase.save_hdf(celldbWithVariability, outputFilename)

# -- Useful for debugging --    
# for k,v in celldbWithTuning.iloc[46].items(): print(f'{k}:\t {v}')
