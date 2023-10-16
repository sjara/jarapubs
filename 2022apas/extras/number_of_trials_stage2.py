"""
Find number of trials during stage 2, so we can match cohorts.

Run as follows (including cohort):
run number_of_trials_stage2.py 3
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
#import studyutils
from importlib import reload
reload(studyparams)

if len(sys.argv)>1:
    COHORT = int(sys.argv[1])
else:
    raise ValueError('You need to specify which cohort to load.')

if COHORT==2:
    '''
    subjectsToInclude = ['pamo009', 'pamo010', 'pamo011', 'pamo012', 'pamo013', 'pamo014',
                         'pamo015', 'pamo017', 'pamo021', 'pamo022', 'pamo024', 'pamo026']
    miceActiveOnly = list(set(studyparams.MICE_ACTIVE_ONLY) & set(subjectsToInclude))
    miceActivePassive = list(set(studyparams.MICE_ACTIVE_PASSIVE) & set(subjectsToInclude))
    #subjects = studyparams.MICE_ALL
    #subjects = miceActiveOnly
    subjects = miceActivePassive
    '''
    subjects = studyparams.MICE_ALL_COH2
    sessionRange = studyparams.SESSIONS_COH2_STAGE2
    excludeSessions = studyparams.EXCLUDE_SESSIONS
elif COHORT==3:
    subjects = studyparams.MICE_ALL_COH3
    #subjects = studyparams.MICE_ACTIVE_ONLY_COH3
    #subjects = studyparams.MICE_PASSIVE_BEFORE_ACTIVE
    sessionRange = studyparams.SESSIONS_COH3_STAGE2
    excludeSessions = studyparams.EXCLUDE_SESSIONS
elif COHORT==4:
    subjects = studyparams.MICE_ALL_COH4
    sessionRange = studyparams.SESSIONS_COH4_STAGE2
    excludeSessions = studyparams.EXCLUDE_SESSIONS

#subjects.remove('pamo025')
#subjects = ['pamo009']
paradigm = '2afc'

#varlist = ['outcomeMode', 'antibiasMode', 'automationMode', 'delayToTarget']
varlist = ['outcomeMode']

sessions = behavioranalysis.sessions_in_range(sessionRange)
dframe = pd.DataFrame(index=sessions)
for indsub, subject in enumerate(subjects):
    toExclude = excludeSessions[subject] if (subject in excludeSessions) else []
    sessions = behavioranalysis.sessions_in_range(sessionRange, toExclude)

    print(f'Loading {subject} {sessionRange}  exclude:{toExclude}')
    
    dflist = []
    for session in sessions:
        behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
        bdata = loadbehavior.BehaviorData(behavFile, varlist)

        #rig = int(bdata.session['hostname'][7:]) # Ignore 'jararig' string
        #delay = bdata['delayToTarget'][-1]
        nTrials = len(bdata['outcomeMode'])
        
        dflist.append({'session':session, subject:nTrials})
        
    sdframe = pd.DataFrame(dflist)
    sdframe = sdframe.set_index('session')

    dframe = dframe.join(sdframe)
    

nTotalEachMouse = dframe.sum(axis=0)

plt.clf()
gsMain = gridspec.GridSpec(1, 6)
gsMain.update(left=0.075, right=0.98, top=0.9, bottom=0.2, wspace=0.5, hspace=0.1)
ax0 = plt.subplot(gsMain[0, 0:-1])
plt.plot(dframe,'o-', color='0.75')
#plt.set_axis
plt.plot(dframe.median(axis=1),'o-', color='C0', lw=4)
plt.ylim([0, 700])
plt.xticks(rotation=80)
plt.title(f'Trials each session (Cohort {COHORT})')

ax1 = plt.subplot(gsMain[0, -1])
plt.plot(np.tile(0,len(nTotalEachMouse)), nTotalEachMouse,'o', color='0.75')
plt.plot([-0.5,0.5], np.tile(nTotalEachMouse.median(),2), 'C0', lw=3)
plt.plot([-0.5,0.5], np.tile(nTotalEachMouse.mean(),2), '0.75', lw=2)
plt.xlim([-1,1])
plt.ylim([1000,5000])
plt.title('Total trials')
plt.show()

sys.exit()


