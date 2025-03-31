"""
Extract running traces from videos and estimate running state on each trial.

You can either specify the subject to analyze, as in:
  python estimate_running.py test000
Or run it without an argument to analyze all subjects:
  python estimate_running.py

If a session has a video but doesn't have a sync light, the onsets are guessed (as periodic).
If the number of sync lights is smaller than the number of trials, NaN is used on those trials.

"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from jaratoolbox import facemapanalysis
from jaratoolbox import celldatabase
from jaratoolbox import loadbehavior
from jaratoolbox import settings
from importlib import reload
import studyparams
reload(facemapanalysis)
reload(studyparams)

SAVE = 1
PLOT = 0

if len(sys.argv)==2:
    subjects = [sys.argv[1]]
else:
    subjects = studyparams.SUBJECTS

PREPOSTSYNC = studyparams.PREPOST_SYNC_LIGHT_START_DATE
outputDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)

sessionTypes = {'PureTones': {'video':'01', 'behavior':'a', 'paradigm':'am_tuning_curve'}}
sessionTypes.update(studyparams.ODDBALL_SESSION_TYPES)

for subject in subjects:
    procDir = os.path.join(settings.VIDEO_PATH, f'{subject}_processed')
    behavDir = os.path.join(settings.BEHAVIOR_PATH, f'{subject}/')

    inforecFile = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')
    inforec = celldatabase.read_inforec(inforecFile)
    sessionDates = [experiment.date for experiment in inforec.experiments]

    runningEachSession = {}  # To store the running value for each trial in each session.

    #for oneDate in [sessionDates[0]]:  #sessionDates: #[sessionDates[-1]]:  #sessionDates: #
    for oneDate in sessionDates:
        for reagent in studyparams.REAGENTS:
            for sessionType, sessionInfo in sessionTypes.items():
                # -- Load the traces from the video (estimated by FaceMap) --
                videoSuffix = sessionInfo['video']
                paradigm = sessionInfo['paradigm']
                behavSuffix = sessionInfo['behavior']

                dateNoDash = oneDate.replace('-','')
                procFilename = f'{subject}_{paradigm}_{dateNoDash}_{videoSuffix}{reagent}_proc.npy'
                procFilepath = os.path.join(procDir, procFilename)
                proc = np.load(procFilepath, allow_pickle=True).item()
                pixchange = proc['pixelchange'][0]

                behavSession = dateNoDash+behavSuffix+reagent
                behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, behavSession)
                behavFilenameOnly = os.path.basename(behavFile)
                bdata = loadbehavior.BehaviorData(behavFile)
                nTrialsBehav = len(bdata['isiMean'])

                syncLight = None
                if proc['blink']:
                    syncLight = proc['blink'][0]
                    if (syncLight.max()-syncLight.min())/syncLight.max() > 0.2:
                        prepost = (oneDate >= PREPOSTSYNC)
                        syncLightOnset = facemapanalysis.find_sync_light_onsets(syncLight,
                                                                                fixmissing=True,
                                                                                prepost=prepost)
                    else:
                       print(f'WARNING! {procFilename} does not seem to have a valid sync light.')
                       syncLight = None

                if syncLight is None:
                    #print(f'No sync light on {procFilename}. Will return array with NaN.')
                    #runningEachSession[behavFilenameOnly] = np.full(nTrialsBehav, np.nan)
                    #continue
                    print(f'No sync light on {procFilename}. We will guess the sync onsets.')
                    nframes = len(pixchange)
                    trialStartTime = bdata.events['eventTime'][bdata.events['eventCode']==-1]
                    trialStartTime = trialStartTime[1:nTrialsBehav+1]
                    syncLightOnset = facemapanalysis.guess_sync_light_onsets(nframes, trialStartTime)


                # -- Account for sync-light happening after stim offset for PureTones session --
                if sessionType=='PureTones':
                    delayToSync = bdata['delayToSyncLight'][-1]
                    if bdata['syncLightMode'][-1] == bdata.labels['syncLightMode']['from_stim_offset']:
                        delayToSync += bdata['stimDur'][-1]
                    syncPeriodInFrames = np.median(np.diff(np.flatnonzero(syncLightOnset)))
                    trialStartTime = bdata.events['eventTime'][bdata.events['eventCode']==-1]
                    trialPeriodInSec = np.median(np.diff(trialStartTime))
                    estimatedFrameRate = syncPeriodInFrames/trialPeriodInSec
                    delayToSyncInFrames = delayToSync * estimatedFrameRate
                    syncLightOnset = np.roll(syncLightOnset, -int(np.round(delayToSyncInFrames)))
                    print(f'Estimated frame rate: {estimatedFrameRate} frames/sec')
                    
                # -- Estimate running on each trial --
                runningEachTrial, runningTraceSmooth = \
                    facemapanalysis.estimate_running_each_trial(pixchange, syncLightOnset,
                                                                smoothsize=10, presamples=4,
                                                                showfig=PLOT)
                # -- Check for 1 extra sync light trial (this is normal) and remove it if present --
                if len(runningEachTrial) == nTrialsBehav+1:
                    runningEachTrial = runningEachTrial[:-1]

                # -- Check that the number of sync lights matches the number of trials --
                trialDiff = nTrialsBehav - len(runningEachTrial)
                errStr = 'WARNING!' if trialDiff else ''
                print(f'Processed {procFilename}: nSync={len(runningEachTrial)}, nTrials={nTrialsBehav}\t{errStr}')

                # -- If there are missing sync lights, fill in the rest assuming same as the last one --
                if trialDiff>0:
                    runningEachTrial = np.concatenate((runningEachTrial, np.tile(np.nan, trialDiff)))
                    print(f'    Fixed: nSync={len(runningEachTrial)}, nTrials={nTrialsBehav}')
                elif trialDiff<0:
                    plt.clf()
                    plt.plot(syncLight)
                    plt.plot(syncLightOnset*np.max(syncLight), '.')
                    plt.title(procFilename)
                    plt.show()
                    plt.xlim([-10, 300])
                    raise ValueError('There are more sync lights than trials. This should not happen.')

                # -- Store running array in dictionary --
                runningEachSession[behavFilenameOnly] = runningEachTrial

                if PLOT:
                    plt.title(behavSession)
                    plt.draw()
                    plt.show()
                    #plt.pause(1)
                    plt.waitforbuttonpress()
                    #sys.exit()
                    #if behavSession == '20230322adoi':
                    #    sys.exit()
                
    outputFilename = os.path.join(outputDir, f'running_{subject}.npz')
    if SAVE:
        np.savez(outputFilename, **runningEachSession)
        print(f'Saved {outputFilename}')

