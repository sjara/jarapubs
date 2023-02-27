"""
Save psychometric peformance for each animal.

Run as:
run -t generate_psychometric_each_mouse.py late 3

"""

import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

if len(sys.argv)>2:
    PERIOD = sys.argv[1]
    COHORT = int(sys.argv[2])
else:
    raise ValueError('You need to specify which period and cohort to process.')

SAVE_RESULTS = 1

if PERIOD=='early':
    dayRange = [0, 4]
if PERIOD=='late':
    #dayRange = [9, 12]
    dayRange = [20, 24]
    #dayRange = [17, 20]

FIGNAME = 'psychometrics_comparison'
figDataFile = f'psychometric_{PERIOD}_coh{COHORT}.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()
figDataFullPath = os.path.join(figDataDir, figDataFile)
scriptFullPath = os.path.realpath(__file__)

varlist = ['outcome', 'valid', 'choice', 'targetFMslope', 'startFreq', 'endFreq', 'targetDuration']
paradigm = '2afc'
subjects = getattr(studyparams, f'MICE_ALL_COH{COHORT}')

nPsySteps = 6  # HARDCODED
avgCorrect = np.empty(len(subjects))
fractionLeft = np.empty((len(subjects), nPsySteps))
ciFractionLeft = np.empty((len(subjects), 2, nPsySteps))
psyCurveParams = np.empty((len(subjects), 4))

for indsub, subject in enumerate(subjects):
    print(f'{subject}')

    sessions = studyutils.get_sessions(subject, stage=4)[slice(*dayRange)]
    if len(sessions)==0:
        print(f'No data available for {subject}. It will be ignored.')
        subjects.remove(subject)
        continue
    
    bdata = behavioranalysis.load_many_sessions(subject, sessions, varlist=varlist)

    correct = bdata['outcome']==bdata.labels['outcome']['correct']
    valid = bdata['valid'].astype(bool)
    leftwardChoice = bdata['choice']==bdata.labels['choice']['left']

    #targetParamValue = bdata['targetFMslope']
    #possibleParamValue = np.unique(targetParamValue)
    #nParamValues = len(possibleParamValue)
    targetDuration = bdata['targetDuration'][-1]
    targetParamValue = np.log2(bdata['endFreq']/bdata['startFreq'])/targetDuration

    (possibleValues, fractionHitsEachValue, ciHitsEachValue, nTrialsEachValue, nHitsEachValue)=\
        behavioranalysis.calculate_psychometric(leftwardChoice, targetParamValue, valid)

    # -- Fit sigmoidal curve --
    par0 = [0, 0.5, 0, 0]
    bounds = [[-np.inf, 0.08, 0, 0], [np.inf, np.inf, 0.5, 0.5]]
    curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues,
                                                  fractionHitsEachValue, p0=par0, bounds=bounds)

    avgCorrect[indsub] = correct.sum()/valid.sum()
    fractionLeft[indsub] = fractionHitsEachValue
    ciFractionLeft[indsub] = ciHitsEachValue
    psyCurveParams[indsub] = curveParams

    if 0:
        plt.cla()
        gsMain = gridspec.GridSpec(1, 1)
        gsMain.update(left=0.15, right=0.95, top=0.9, bottom=0.15)
        ax0 = plt.subplot(gsMain[0, 0])

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
        plt.xlim([-6, 6]) #plt.xlim([-1.2, 1.2])
        plt.ylabel('Leftward choice (%)', fontsize=12)
        plt.xlabel('FM slope (oct/s)', fontsize=12)
        titleStr = f'Psychometric performance, one session'
        plt.title(titleStr, fontsize=12, fontweight='bold')
        plt.text(1, 10, f'{subject}\n{sessions}', color='0.9', ha='right')
        extraplots.set_ticks_fontsize(plt.gca(), 12)
        extraplots.boxoff(ax0)
        #plt.waitforbuttonpress()
        plt.pause(0.5)
        plt.show()

if SAVE_RESULTS:
    np.savez(figDataFullPath, subjects=subjects, dayRange=dayRange, avgCorrect=avgCorrect,
             possibleFMslope=possibleValues, fractionLeft=fractionLeft,
             ciFractionLeft=ciFractionLeft, psyCurveParams=psyCurveParams)
    print(f'Saved {figDataFullPath}')
