"""
Find number of trials per stage.

Run script indicating cohort to process:
run generate_trials_per_stage.py 4

Note that you first need to run generate_sessions_each_mouse.py for that cohort.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from jaratoolbox import settings
from jaratoolbox import loadbehavior
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import studyparams
import studyutils
import figparams

if len(sys.argv)>1:
    COHORT = int(sys.argv[1])
else:
    raise ValueError('You need to specify which cohort to process.')

SAVE_RESULTS = 0

FIGNAME = 'trials_per_stage'
figDataFile = f'trials_per_stage_coh{COHORT}.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()
figDataFullPath = os.path.join(figDataDir, figDataFile)
scriptFullPath = os.path.realpath(__file__)

subjects = getattr(studyparams, f'MICE_ALL_COH{COHORT}')
#subjects = getattr(studyparams, f'MICE_ACTIVE_ONLY_COH2')
#subjects = getattr(studyparams, f'MICE_ACTIVE_ONLY_COH3')
#subjects = getattr(studyparams, f'MICE_ACTIVE_PASSIVE') # Coh2
#subjects = getattr(studyparams, f'MICE_PASSIVE_THEN_ACTIVE') # Coh3
### BADsubjects = [['pamo009', 'pamo011', 'pamo012', 'pamo013', 'pamo022', 'pamo026',
###             'pamo010', 'pamo014', 'pamo015', 'pamo017', 'pamo021', 'pamo024']] # Fast Coh2
#subjects = ['pamo028', 'pamo032', 'pamo039', 'pamo040', 'pamo041'] # Fast Coh3

#subjects.remove('pamo052')

dflist = []  # To store data for dataframe

for indsub, subject in enumerate(subjects):
    print(f'{indsub}  {subject}')
    dframe = studyutils.load_sessions_dataframe(subject)

    # Fixing bug: when the paradigm stops automatically, it does not store th last trial.
    # Because of this, we need to add one trial (only) to sessions with 499 trials.
    dframe.ntrials += (dframe.ntrials == 499).astype(int)

    # Insert estimates from excluded sessions
    if subject in studyparams.EXCLUDE_SESSIONS:
        for onedate in studyparams.EXCLUDE_SESSIONS[subject]:
            dayBefore, thisDay, dayAfter = studyutils.days_around(onedate, outformat='%Y%m%da')
            #print( dayBefore, thisDay, dayAfter )
            entryBefore = dframe[dframe.session==dayBefore]
            entryAfter = dframe[dframe.session==dayAfter]
            newEntry = entryBefore.copy()
            newEntry.session = thisDay
            newEntry.ntrials = (int(entryBefore.ntrials) + int(entryAfter.ntrials))//2
            dframe = dframe.append(newEntry, ignore_index=True)
    
    meanPerStage = dframe.groupby('stage').mean()['ntrials'] # Mean trials per session each stage
    sumPerStage = dframe.groupby('stage').sum()['ntrials'] # Total trials for each stage

    totalPreS3 = sumPerStage[0:3].sum() ; plotLabel='Total S0-S2'
    totalPreS3 = sumPerStage[0:2].sum(); print('Only S0-S1'); plotLabel='Total S0-S1'
    totalPreS3 = sumPerStage[2].sum(); print('Only S2'); plotLabel='Total S2'

    if len(meanPerStage)>3:
        dfDict = {'totalPreS3':totalPreS3, 'totalS3':sumPerStage[3], 'meanS3':meanPerStage[3]}
    else:
        dfDict = {'totalPreS3':totalPreS3, 'totalS3':0, 'meanS3':0}        
    dflist.append(dfDict)
    
dfTrialsPerStage = pd.DataFrame(dflist, index=subjects)

if 1:
    plt.clf()
    plt.subplot(1,3,1)
    plt.plot(np.tile(0,len(subjects)), dfTrialsPerStage.totalPreS3, 'o', mfc='none')
    print(f'Median S0-S2: {dfTrialsPerStage.totalPreS3.median()}')
    #plt.ylim([0, 4000])
    plt.ylim([0, 8000])
    plt.ylabel(plotLabel)  #'Total S0-S2'
    plt.subplot(1,3,2)
    plt.plot(np.tile(0,len(subjects)), dfTrialsPerStage.totalS3, 'o', mfc='none')
    plt.ylim([10000, 18000])
    plt.ylabel('Total S3')
    plt.subplot(1,3,3)
    plt.plot(np.tile(0,len(subjects)), dfTrialsPerStage.meanS3, 'o', mfc='none')
    plt.ylim([400, 520])
    plt.ylabel('Trials per session S3')
    plt.show()

if SAVE_RESULTS:
    dfTrialsPerStage.to_csv(figDataFullPath)
    print(f'Saved {figDataFullPath}')


