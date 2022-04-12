"""
Find sessions for each mouse on each stage.
"""

import sys
import os
import pandas as pd
import studyparams
import jaratoolbox
from jaratoolbox import loadbehavior
from jaratoolbox import behavioranalysis
from importlib import reload
reload(studyparams)

outputDir = '/tmp/sessions/'
#outputDir = os.path.join(jaratoolbox.settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'sessions')

sessionRange = studyparams.SESSIONS_RANGE
excludeSessions = studyparams.EXCLUDE_SESSIONS

subjects = studyparams.MICE_ALL
paradigm = '2afc'

varlist = ['outcomeMode', 'antibiasMode', 'psycurveMode']

for indsub, subject in enumerate(subjects):
    toExclude = excludeSessions[subject] if (subject in excludeSessions) else []
    sessions = behavioranalysis.sessions_in_range(sessionRange, toExclude)

    print(f'Loading {subject} {sessionRange}  exclude:{toExclude}')
    
    dflist = []
    for session in sessions:
        behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
        bdata = loadbehavior.BehaviorData(behavFile, varlist)
        if bdata['outcomeMode'][-1] == bdata.labels['outcomeMode']['sides_direct']:
            stage = 0
        elif bdata['outcomeMode'][-1] == bdata.labels['outcomeMode']['direct']:
            stage = 1
        elif bdata['outcomeMode'][-1] == bdata.labels['outcomeMode']['on_next_correct']:
            stage = 2
        elif bdata['outcomeMode'][-1] == bdata.labels['outcomeMode']['only_if_correct']:
            if bdata['psycurveMode'][-1] == bdata.labels['psycurveMode']['off']:
                stage = 3
            elif bdata['psycurveMode'][-1] == bdata.labels['psycurveMode']['uniform']:
                stage = 4
            elif bdata['psycurveMode'][-1] == bdata.labels['psycurveMode']['controls']:
                stage = 5
        else:
            stage = None

        if bdata['antibiasMode'][-1] == bdata.labels['antibiasMode']['repeat_mistake']:
            antibias = 1
        elif bdata['antibiasMode'][-1] == bdata.labels['antibiasMode']['off']:
            antibias = 0
        else:
            antibias = None

        dflist.append({'session':session, 'stage':stage, 'antibias':antibias})

    dframe = pd.DataFrame(dflist)

    outputFile = os.path.join(outputDir, f'{subject}.csv')
    dframe.to_csv(outputFile)
    print(f'Saved {outputFile}')

