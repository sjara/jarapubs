"""
Save data for learning curve stage 3.

Run this script with an argument specifying which cohort to evaluate:
run -t generate_learning_curve_stage3 4
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

SAVE_RESULTS = 1

FIGNAME = 'learning_curve_stage3'
figDataFile = f'fraction_correct_stage3_coh{COHORT}.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()
figDataFullPath = os.path.join(figDataDir, figDataFile)
scriptFullPath = os.path.realpath(__file__)

varlist = ['outcome', 'valid', 'choice']
paradigm = '2afc'
subjects = getattr(studyparams, f'MICE_ALL_COH{COHORT}')

dflist = []  # To store data for dataframe
plt.clf()
for indsub, subject in enumerate(subjects):
    print(f'{indsub}  {subject}')
    sessions = studyutils.get_sessions(subject, stage=3)
    bdata = behavioranalysis.load_many_sessions(subject, sessions, paradigm, varlist=varlist)

    nSessions = len(sessions)
    assert nSessions == (bdata['sessionID'][-1] + 1)
    correct = bdata['outcome']==bdata.labels['outcome']['correct']
    valid = bdata['valid'].astype(bool) & (bdata['choice']!=bdata.labels['choice']['none'])

    correctEachSession = np.empty(nSessions)
    validEachSession = np.empty(nSessions)
    for sessionID in range(nSessions):
        trialsThisSession = bdata['sessionID']==sessionID
        correctEachSession[sessionID] = np.sum(correct[trialsThisSession])
        validEachSession[sessionID] = np.sum(valid[trialsThisSession])

    fractionCorrect = correctEachSession/validEachSession

    dflist.append(dict(zip(sessions, fractionCorrect)))
    
    if 1:
        plt.plot(fractionCorrect,'.-', lw=2, ms=16)
        plt.ylim([0.40, 1.0])
        plt.ylabel('Fraction correct')
        plt.xlabel('Session')
        plt.grid(True)
        plt.title(subject)
        plt.pause(0.01)
        plt.show()

dfFractionCorrect = pd.DataFrame(dflist, index=subjects)

if SAVE_RESULTS:
    dfFractionCorrect.to_csv(figDataFullPath)
    print(f'Saved {figDataFullPath}')


'''
plt.clf()
colorAO = figparams.colors['activeOnly']
colorAP = figparams.colors['activePassive']

dfAO = dframe.loc[studyparams.MICE_ACTIVE_ONLY]
dfAP = dframe.loc[studyparams.MICE_ACTIVE_PASSIVE]
#plt.plot(dframe.median(),'o-', color='k', lw=3); ylim([0.5, 0.9])
plt.plot(dfAO.median(),'o-', color=colorAO, lw=3);
plt.plot(dfAP.median(),'o-', color=colorAP, lw=3);
plt.ylim([0.5, 0.9])
plt.xticks(rotation=45)
plt.show()
'''

'''
dfFractionCorrect = pd.read_csv('/data/figuresdata/2022apas/learning_curve_stage3/fraction_correct_stage3_coh3.csv', index_col=0)
for subject, fractionCorrect in dfFractionCorrect.iterrows():
    plt.cla()
    plt.plot(fractionCorrect,'.-', lw=2, ms=16)
    plt.ylim([0.40, 1.0])
    plt.ylabel('Fraction correct')
    plt.xlabel('Session')
    plt.grid(True)
    plt.title(subject)
    plt.waitforbuttonpress()
    plt.show()
'''

