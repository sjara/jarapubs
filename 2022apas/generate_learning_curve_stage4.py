"""
Save data for learning curve stage 4.
"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.optimize
from jaratoolbox import settings
from jaratoolbox import loadbehavior
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import studyparams
import studyutils
import figparams
from importlib import reload
reload(studyparams)

if len(sys.argv)>1:
    COHORT = int(sys.argv[1])
else:
    raise ValueError('You need to specify which cohort to process.')

SAVE_RESULTS = 0

FIGNAME = 'learning_curve_stage4'
figDataFile = f'measurements_stage4_coh{COHORT}.csv'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()
figDataFullPath = os.path.join(figDataDir, figDataFile)
scriptFullPath = os.path.realpath(__file__)

varlist = ['outcome', 'valid', 'choice']
paradigm = '2afc'
subjects = getattr(studyparams, f'MICE_ALL_COH{COHORT}')

#subjects = ['pamo009']
#subjects = ['pamo026']
#subjects = ['pamo025', 'pamo026'] #
#subjects = studyparams.MICE_ALL

measuresLabels = ['bias', 'invSlope', 'lapseHigh', 'lapseLow', 'nValid', 'nCorrect']
nMeasures = len(measuresLabels)
dflists = [[] for ind in range(nMeasures)]
perflist = []

plt.clf()
for indsub, subject in enumerate(subjects):
    print(f'{indsub}  {subject}')
    sessions = studyutils.get_sessions(subject, stage=4)
    ###sessions = sessions[:1] ######## TESTING ########
    nSessions = len(sessions)
    
    curveParamsThisSubject = np.empty((nSessions, 4))
    correctEachSession = np.empty(nSessions)
    validEachSession = np.empty(nSessions)
    for indsession, session in enumerate(sessions):
        behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
        bdata = loadbehavior.BehaviorData(behavFile)

        correct = bdata['outcome']==bdata.labels['outcome']['correct']
        valid = bdata['valid'].astype(bool)
        rightwardChoice = bdata['choice']==bdata.labels['choice']['right']

        correctEachSession[indsession] = np.sum(correct)
        validEachSession[indsession] = np.sum(valid)
        
        targetParamValue = bdata['targetFMslope']
        possibleParamValue = np.unique(targetParamValue)
        nParamValues = len(possibleParamValue)

        (possibleValues, fractionHitsEachValue, ciHitsEachValue, nTrialsEachValue, nHitsEachValue)=\
            behavioranalysis.calculate_psychometric(rightwardChoice, targetParamValue, valid)

        # -- Fit sigmoidal curve --
        par0 = [0, -0.5, 0, 0]
        bounds = [[-np.inf, -np.inf, 0, 0], [np.inf, -0.08, 0.5, 0.5]]
        curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues,
                                                      fractionHitsEachValue, p0=par0, bounds=bounds)

        print(f'{subject} {session}: {curveParams}')
        curveParamsThisSubject[indsession,:] = curveParams 

        if 1:
            plt.cla()
            xPad = 0.2 * (possibleValues[-1] - possibleValues[0])
            fitxval = np.linspace(possibleValues[0]-xPad, possibleValues[-1]+xPad, 40)
            fityval = extrastats.psychfun(fitxval, *curveParams)
            xTicks = np.arange(-1, 1.5, 0.5)
            fontSizeLabels = 12
            hfit = plt.plot(fitxval, 100*fityval, '-', linewidth=2, color='k')
            (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValues,
                                                                        fractionHitsEachValue,
                                                                        ciHitsEachValue, xTicks=None,
                                                                        xscale='linear')
            pline.set_visible(False)
            plt.xlim([-1.2, 1.2])
            plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
            plt.xlabel('FM slope (A.U.)', fontsize=fontSizeLabels)
            titleStr = f'{subject}: {session}'
            plt.title(titleStr, fontsize=fontSizeLabels, fontweight='bold')
            plt.grid(True, axis='y', color='0.9')
            plt.show()
            plt.pause(0.01)
            plt.waitforbuttonpress()

    measuresThisSuject = np.hstack((curveParamsThisSubject,
                                    validEachSession[:,np.newaxis],
                                    correctEachSession[:,np.newaxis]))
    for indp, param in enumerate(measuresLabels):
        #dflists[indp].append(dict(zip(sessions, curveParamsThisSubject[:,indp])))
        dflists[indp].append(dict(zip(sessions, measuresThisSuject[:,indp])))
    #fractionCorrect = correctEachSession/validEachSession
    #perflist.append(dict(zip(sessions,fractionCorrect)))

measureDict = {}
for indp, measure in enumerate(measuresLabels):
    dfThisMeasure = pd.DataFrame(dflists[indp], index=subjects)
    measureDict[measure] = dfThisMeasure

dfMeasures = pd.concat(measureDict)
if SAVE_RESULTS:
    dfMeasures.to_csv(figDataFullPath)
    print(f'Saved {figDataFullPath}')

#dfFractionCorrect = pd.DataFrame(perflist, index=subjects)

# -- To access multi-index --
#  dfCurveParams.xs('invSlope')
#  dfCurveParams.xs('pamo025', level=1)
