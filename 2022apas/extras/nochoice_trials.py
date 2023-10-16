"""
Estimate number of trials with no choice.
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
sys.path.append('..')
import studyparams
import studyutils
import figparams


varlist = ['outcome', 'valid', 'choice']
paradigm = '2afc'
subjects = (getattr(studyparams, 'MICE_ALL_COH2') +
            getattr(studyparams, 'MICE_ALL_COH3'))

nTrials = np.empty(len(subjects), dtype=int)
nNoChoice = np.empty(len(subjects), dtype=int)

for indsub, subject in enumerate(subjects):
    print(f'{indsub}  {subject}')
    sessions = studyutils.get_sessions(subject, stage=3)
    bdata = behavioranalysis.load_many_sessions(subject, sessions, paradigm, varlist=varlist)

    nSessions = len(sessions)
    assert nSessions == (bdata['sessionID'][-1] + 1)
    noChoice = bdata['outcome']==bdata.labels['outcome']['nochoice']
    #noChoice2 = bdata['choice']==bdata.labels['choice']['none']

    nTrials[indsub] = len(noChoice)
    nNoChoice[indsub] = sum(noChoice)


fractionNoChoice = nNoChoice/nTrials
meanNoChoice = np.mean(fractionNoChoice)
stdNoChoice = np.std(fractionNoChoice)

print(f'No choice: {meanNoChoice:0.2%} +/- {stdNoChoice:0.2%}')
