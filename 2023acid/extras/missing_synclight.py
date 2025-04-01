"""
Test what to do for session with missing sync light.
"""

import os
import sys
sys.path.append('..')
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

subject = 'acid006'

procDir = os.path.join(settings.VIDEO_PATH, f'{subject}_processed')
behavDir = os.path.join(settings.BEHAVIOR_PATH, f'{subject}/')

inforecFile = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')
inforec = celldatabase.read_inforec(inforecFile)

#sessionDates = [experiment.date for experiment in inforec.experiments]

reagent = studyparams.REAGENTS[1]  # 'pre', 'saline', 'doi'
#sessionType = {'PureTones': {'video':'01', 'behavior':'a', 'paradigm':'am_tuning_curve'}}
sessionInfo = {'video':'01', 'behavior':'a', 'paradigm':'am_tuning_curve'}
oneDate = '2023-03-22'

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

if proc['blink']:
    syncLight = proc['blink'][0]

trialStartTime = bdata.events['eventTime'][bdata.events['eventCode']==-1]
nTrials = len(bdata['stimDur'])
trialStartTime = trialStartTime[1:nTrials+1]
nFrames = len(pixchange)

syncLightOnset = facemapanalysis.guess_sync_light_onsets(nFrames, trialStartTime, framerate=29)

plt.clf()
plt.plot(syncLightOnset)
plt.show()

