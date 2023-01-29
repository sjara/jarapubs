"""
Test performance of pamo mice under the control FM stimuli.
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


paradigm = '2afc'
subjects = studyparams.MICE_ALL_COH2

#subjects = [subjects[0]]

sessions = ['20220406a']


plt.clf()
for indsub, subject in enumerate(subjects):
    print(f'{indsub}  {subject}')

    #sessions = studyutils.get_sessions(subject, stage=4)
    nSessions = len(sessions)
    
    #curveParamsThisSubject = np.empty((nSessions, 4))
    
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
        startFreq = bdata['startFreq']
        possibleStartFreq = np.unique(startFreq)
        endFreq = bdata['endFreq']
        possibleEndFreq = np.unique(endFreq)

        '''
        hiSubset = ( (targetParamValue==possibleParamValue[0]) |
                     ((targetParamValue==possibleParamValue[1]) & (startFreq==possibleStartFreq[0])) |
                     ((targetParamValue==possibleParamValue[2]) & (startFreq==possibleStartFreq[1])) |
                     (targetParamValue==possibleParamValue[3]) )
        '''
        extSubset = (startFreq==possibleStartFreq[0]) | (startFreq==possibleStartFreq[2])
        #extSubset = (endFreq==possibleEndFreq[1])
        midSubset = (startFreq==possibleStartFreq[1])
        
        (possibleValues1, fractionHitsEachValue1, ciHitsEachValue1, nTrialsEachValue1, nHitsEachValue1)=\
            behavioranalysis.calculate_psychometric(rightwardChoice, targetParamValue, valid & extSubset)
        (possibleValues2, fractionHitsEachValue2, ciHitsEachValue2, nTrialsEachValue2, nHitsEachValue2)=\
            behavioranalysis.calculate_psychometric(rightwardChoice, targetParamValue, valid & midSubset)

        # -- Fit sigmoidal curve --
        '''
        par0 = [0, -0.5, 0, 0]
        bounds = [[-np.inf, -np.inf, 0, 0], [np.inf, -0.08, 0.5, 0.5]]
        curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues,
                                                      fractionHitsEachValue, p0=par0, bounds=bounds)

        print(f'{subject} {session}: {curveParams}')
        curveParamsThisSubject[indsession,:] = curveParams 
        '''
        
        if 1:
            plt.cla()
            #xPad = 0.2 * (possibleValues[-1] - possibleValues[0])
            #fitxval = np.linspace(possibleValues[0]-xPad, possibleValues[-1]+xPad, 40)
            #fityval = extrastats.psychfun(fitxval, *curveParams)
            xTicks = np.arange(-1, 1.5, 0.5)
            fontSizeLabels = 12
            #hfit = plt.plot(fitxval, 100*fityval, '-', linewidth=2, color='k')
            (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValues1,
                                                                        fractionHitsEachValue1,
                                                                        ciHitsEachValue1, xTicks=None,
                                                                        xscale='linear')
            (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValues2,
                                                                        fractionHitsEachValue2,
                                                                        ciHitsEachValue2, xTicks=None,
                                                                        xscale='linear')
            #pline.set_visible(False)
            plt.xlim([-1.2, 1.2])
            plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
            plt.xlabel('FM slope (A.U.)', fontsize=fontSizeLabels)
            titleStr = f'{subject}: {session}'
            plt.title(titleStr, fontsize=fontSizeLabels, fontweight='bold')
            plt.grid(True, axis='y', color='0.9')
            plt.show()
            plt.pause(0.01)
            plt.waitforbuttonpress()

    '''
    measuresThisSuject = np.hstack((curveParamsThisSubject,
                                    validEachSession[:,np.newaxis],
                                    correctEachSession[:,np.newaxis]))
    for indp, param in enumerate(measuresLabels):
        #dflists[indp].append(dict(zip(sessions, curveParamsThisSubject[:,indp])))
        dflists[indp].append(dict(zip(sessions, measuresThisSuject[:,indp])))
    #fractionCorrect = correctEachSession/validEachSession
    #perflist.append(dict(zip(sessions,fractionCorrect)))
    '''

'''
measureDict = {}
for indp, measure in enumerate(measuresLabels):
    dfThisMeasure = pd.DataFrame(dflists[indp], index=subjects)
    measureDict[measure] = dfThisMeasure

dfMeasures = pd.concat(measureDict)

#dfMeasures.to_csv(figDataFullPath)
#print(f'Saved {figDataFullPath}')
'''
