"""
Check criteria for each stage. 
For example in S3 if animals are biased on the last session to set anti-bias on next session.

USAGE:
run extra_check_stage_criteria.py 20220628
"""

import sys
import os
import numpy as np
import pandas as pd
import jaratoolbox
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import loadbehavior
from jaratoolbox import behavioranalysis
sys.path.append('..')
import studyparams
import studyutils
from importlib import reload
reload(studyparams)
reload(studyutils)


if len(sys.argv)==2:
    session = sys.argv[1]+'a'
else:
    session = '20220628a'

subjects = studyparams.MICE_ALL_COH4
#subjects = studyparams.MICE_ALL_COH3
#subjects = studyparams.MICE_ALL_COH2
#subjects.remove('pamo052')
#subjects = ['test000']

paradigm = '2afc'
varlist = ['antibiasMode', 'outcome', 'valid', 'choice', 'rewardSide']

dflist = []
for indsub, subject in enumerate(subjects):
    if subject in studyparams.EXCLUDE_SESSIONS:
        isodate = studyutils.session_to_isodate(session)
        if isodate in studyparams.EXCLUDE_SESSIONS[subject]:
            print(f'Excluding data for {subject}')
            continue
    
    behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
    bdata = loadbehavior.BehaviorData(behavFile, varlist)

    if bdata.session['hostname'][:7]=='jararig':
        rig = int(bdata.session['hostname'][7:]) # Ignore 'jararig' string
    else:
        rig = 0
    #delay = bdata['delayToTarget'][-1]
    correct = bdata['outcome']==bdata.labels['outcome']['correct']
    afterError = bdata['outcome']==bdata.labels['outcome']['aftererror']
    valid = bdata['valid'].astype(bool) & (bdata['choice']!=bdata.labels['choice']['none'])
    rewardSideRight = bdata['rewardSide']==bdata.labels['rewardSide']['right']
    rewardSideLeft = bdata['rewardSide']==bdata.labels['rewardSide']['left']
    rightChoice = bdata['choice']==bdata.labels['choice']['right']
    leftChoice = bdata['choice']==bdata.labels['choice']['right']
    percentCorrect = 100*correct.sum()/valid.sum()
    percentCorrectRight = 100*np.sum(correct & rewardSideRight)/np.sum(valid & rewardSideRight)
    percentCorrectLeft = 100*np.sum(correct & rewardSideLeft)/np.sum(valid & rewardSideLeft)
    nTrials = len(bdata['outcome'])
    nRewards = correct.sum()
    withAfterError = correct.sum() + afterError.sum()
    antibiasMode = bdata['antibiasMode'][-1]
    nextDayAntibiasLessThan20 = (percentCorrectRight<20)|(percentCorrectLeft<20)
    nextDayAntibiasLessThan30 = (percentCorrectRight<30)|(percentCorrectLeft<30)
    
    dflist.append({'subject':subject,
                   '<20%':nextDayAntibiasLessThan20,
                   '<30%':nextDayAntibiasLessThan30,
                   'ABmode':antibiasMode,
                   'nTrials':nTrials,
                   'nRewards':nRewards,
                   'withAfterE': withAfterError,
                   'percentCorrect':percentCorrect,
                   'percentCorrectLeft':percentCorrectLeft,
                   'percentCorrectRight':percentCorrectRight,
                   'rig':rig})

print(f'--- {session} ---')
dframe = pd.DataFrame(dflist)

print(dframe)
print(f'Average performance across mice: {dframe.percentCorrect.mean()}')
#dframe[dframe.nextDayAntibias]
#dframe[~dframe.nextDayAntibias]

sys.exit()


