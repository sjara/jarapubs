"""
Save data for example mouse psychometric- one for VOT animal, one for FT animal
"""

import os
import sys
import numpy as np
from jaratoolbox import settings
from jaratoolbox import behavioranalysis
from jaratoolbox import loadbehavior
from jaratoolbox import extrastats
from jaratoolbox import extraplots
import scipy.optimize
import matplotlib.pyplot as plt
from statsmodels.stats.proportion import proportion_confint #Used to compute confidence interval for the error bars.
import studyparams

FIGNAME = 'figure_behavior'
figDataFile = 'data_example_mice_psycurves.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

# -- Optional --
if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)
    #print('Please create folder: {}'.format(figDataDir)); sys.exit()

figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

#databaseName = 'THIS_STUDY_database.h5'
#databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
#celldb = celldatabase.load_hdf(databaseFullPath)

# -- Deal with the stuff above this! --
paradigm = '2afc_speech'
tasks = ['VOT', 'FT']
allFeatValues = np.zeros((2,6))
allFractionHitsEachValue = np.zeros((2,6))
allCiHitsEachValue = np.zeros((2,6))
allnTrialsEachValue = np.zeros((2,6))
allnHitsEachValue = np.zeros((2,6))


for indTask, thisTask in enumerate(tasks):

    subject = studyparams.EXAMPLE_MICE[thisTask]
    sessionRange = studyparams.SESSIONS_RANGE[thisTask]
    sessions = behavioranalysis.sessions_in_range(sessionRange)

    bdata = behavioranalysis.load_many_sessions(subject, sessions, paradigm)

    rightChoice = bdata['choice'] == bdata.labels['choice']['right']
    valid=bdata['valid']& (bdata['choice']!=bdata.labels['choice']['none'])
    if thisTask == 'FT':
        targetFrequency = bdata['targetFTpercent']
    elif thisTask == 'VOT':
        targetFrequency = bdata['targetVOTpercent']

    (possibleValues, fractionHitsEachValue, ciHitsEachValue, nTrialsEachValue, nHitsEachValue) = behavioranalysis.calculate_psychometric(rightChoice, targetFrequency, valid)
    curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues,
                                                  fractionHitsEachValue)


    if thisTask == 'FT':
        possibleValuesFt = possibleValues
        fractionHitsEachValueFt = fractionHitsEachValue
        ciHitsEachValueFt = ciHitsEachValue
        nTrialsEachValueFt = nTrialsEachValue
        nHitsEachValueFt = nHitsEachValue
        curveParamsFt = curveParams
        pcovFt = pcov
    elif thisTask == 'VOT':
        possibleValuesVot = possibleValues
        fractionHitsEachValueVot = fractionHitsEachValue
        ciHitsEachValueVot = ciHitsEachValue
        nTrialsEachValueVot = nTrialsEachValue
        nHitsEachValueVot = nHitsEachValue
        curveParamsVot = curveParams
        pcovVot = pcov


if 1:
    xPadFt = 0.2 * (possibleValuesFt[-1] - possibleValuesFt[0])
    fitxvalFt = np.linspace(possibleValuesFt[0]-xPadFt, possibleValuesFt[-1]+xPadFt, 40)
    fityvalFt = extrastats.psychfun(fitxvalFt, *curveParamsFt)
    xTicks = np.arange(-1, 1.5, 0.5)
    fontSizeLabels = 12
    hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
    (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesFt,
                                                                fractionHitsEachValueFt,
                                                                ciHitsEachValueFt, xTicks=None,
                                                                xscale='linear')
    pline.set_visible(False)
    plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
    plt.xlabel('FT slope (A.U.)', fontsize=fontSizeLabels)
    plt.title('FT example mouse', fontsize=fontSizeLabels, fontweight='bold')
    plt.grid(True, axis='y', color='0.9')

    plt.show()

    plt.figure()
    xPadVot = 0.2 * (possibleValuesVot[-1] - possibleValuesVot[0])
    fitxvalVot = np.linspace(possibleValuesVot[0]-xPadVot, possibleValuesVot[-1]+xPadVot, 40)
    fityvalVot = extrastats.psychfun(fitxvalVot, *curveParamsVot)
    fontSizeLabels = 12
    hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
    (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesVot,
                                                                fractionHitsEachValueVot,
                                                                ciHitsEachValueVot, xTicks=None,
                                                                xscale='linear')
    pline.set_visible(False)
    plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
    plt.xlabel('Vot (A.U.)', fontsize=fontSizeLabels)
    plt.title('VOT example mouse', fontsize=fontSizeLabels, fontweight='bold')
    plt.grid(True, axis='y', color='0.9')

    plt.show()

np.savez(figDataFullPath, possibleValuesFt = possibleValuesFt, fractionHitsEachValueFt = fractionHitsEachValueFt, ciHitsEachValueFt = ciHitsEachValueFt, nTrialsEachValueFt = nTrialsEachValueFt, nHitsEachValueFt = nHitsEachValueFt, curveParamsFt = curveParamsFt, pcovFt = pcovFt, possibleValuesVot = possibleValuesVot, fractionHitsEachValueVot = fractionHitsEachValueVot, ciHitsEachValueVot = ciHitsEachValueVot, nTrialsEachValueVot = nTrialsEachValueVot, nHitsEachValueVot = nHitsEachValueVot, curveParamsVot = curveParamsVot, pcovVot = pcovVot)
print('saved to ' f'{figDataFullPath}')
