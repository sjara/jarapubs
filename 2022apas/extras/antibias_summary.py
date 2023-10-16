"""
Calculate how many sessions of antibias were used for each cohort.
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


#varlist = ['outcome', 'valid', 'choice', 'antibiasMode']
varlist = ['antibiasMode']
paradigm = '2afc'
#subjects = (getattr(studyparams, 'MICE_ALL_COH2') +
#            getattr(studyparams, 'MICE_ALL_COH3'))
subjects = getattr(studyparams, 'MICE_ACTIVE_ONLY_COH2')
#subjects = getattr(studyparams, 'MICE_ACTIVE_PASSIVE')
subjects = getattr(studyparams, 'MICE_PASSIVE_THEN_ACTIVE')


#nTrials = np.empty(len(subjects), dtype=int)
#nNoChoice = np.empty(len(subjects), dtype=int)
antibiasMode = np.empty(len(subjects), dtype=int)
nSessions = np.empty(len(subjects), dtype=int)

for indsub, subject in enumerate(subjects):
    #print(f'{indsub}  {subject}')
    sessions = studyutils.get_sessions(subject, stage=3)
    #bdata = behavioranalysis.load_many_sessions(subject, sessions, paradigm, varlist=varlist)

    nAntibias = 0
    for session in sessions:
        behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
        bdata = loadbehavior.BehaviorData(behavFile, varlist=varlist)
        nAntibias += bdata['antibiasMode'][-1]
    
    nSessions[indsub] = len(sessions)
    #assert nSessions == (bdata['sessionID'][-1] + 1)
    #noChoice = bdata['outcome']==bdata.labels['outcome']['nochoice']
    antibiasMode[indsub] = nAntibias
    print(f'{subject}  {nAntibias}/{nSessions[indsub]} = {nAntibias/nSessions[indsub]:0.1%} sessions with antibias mode')
    
    #noChoice2 = bdata['choice']==bdata.labels['choice']['none']
    #nTrials[indsub] = len(noChoice)
    #nNoChoice[indsub] = sum(noChoice)


#fractionNoChoice = nNoChoice/nTrials
#meanNoChoice = np.mean(fractionNoChoice)
#stdNoChoice = np.std(fractionNoChoice)

medianAntibias = np.median(antibiasMode)
print(f'Antibias median: {medianAntibias}')
