"""
Save data for cohort average psychometric curves- one for VOT, one for FT
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
figDataFile = 'data_cohort_average_psycurves.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

# -- Optional --
#if not os.path.exists(figDataDir):
#    os.mkdir(figDataDir)
    #print('Please create folder: {}'.format(figDataDir)); sys.exit()

figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)


paradigm = '2afc_speech'
tasks = ['VOT', 'FT']
allFeatValues = np.zeros((2,6))
allFractionHitsEachValue = np.zeros((2,6))
allCiHitsEachValue = np.zeros((2,6))
allnTrialsEachValue = np.zeros((2,6))
allnHitsEachValue = np.zeros((2,6))
possibleValuesFt = np.empty((np.size(studyparams.BEHAVIOR_MICE['FT']),6))
fractionHitsEachValueFt = np.empty((np.size(studyparams.BEHAVIOR_MICE['FT']),6))
ciHitsEachValueFt = np.empty((np.size(studyparams.BEHAVIOR_MICE['FT']),2,6))
nTrialsEachValueFt = np.empty((np.size(studyparams.BEHAVIOR_MICE['FT']),6))
nHitsEachValueFt = np.empty((np.size(studyparams.BEHAVIOR_MICE['FT']),6))
possibleValuesVot = np.empty((np.size(studyparams.BEHAVIOR_MICE['VOT']),6))
fractionHitsEachValueVot = np.empty((np.size(studyparams.BEHAVIOR_MICE['VOT']),6))
ciHitsEachValueVot = np.empty((np.size(studyparams.BEHAVIOR_MICE['VOT']),2,6))
nTrialsEachValueVot = np.empty((np.size(studyparams.BEHAVIOR_MICE['VOT']),6))
nHitsEachValueVot = np.empty((np.size(studyparams.BEHAVIOR_MICE['VOT']),6))
nFTmice = len(studyparams.BEHAVIOR_MICE['FT'])
nVOTmice = len(studyparams.BEHAVIOR_MICE['VOT'])


for indTask, thisTask in enumerate(tasks):
    cohort = studyparams.BEHAVIOR_MICE[thisTask]
    for indMouse, thisMouse in enumerate(cohort):

        subject = thisMouse
        sessionRange = studyparams.SESSIONS_RANGE[thisMouse]
        sessions = behavioranalysis.sessions_in_range(sessionRange)

        bdata = behavioranalysis.load_many_sessions(subject, sessions, paradigm)
        if indMouse == 0:
            rightChoice = bdata['choice'] == bdata.labels['choice']['right']
            valid = bdata['valid']& (bdata['choice'] != bdata.labels['choice']['none'])
            if thisTask == 'FT':
                targetFrequency = bdata['targetFTpercent']
            elif thisTask == 'VOT':
                targetFrequency = bdata['targetVOTpercent']
        else:
            rightChoicethisMouse = bdata['choice'] == bdata.labels['choice']['right']
            validThisMouse = valid = bdata['valid']& (bdata['choice'] != bdata.labels['choice']['none'])
            if thisTask == 'FT':
                targetFrequencythisMouse = bdata['targetFTpercent']
            elif thisTask == 'VOT':
                targetFrequencythisMouse = bdata['targetVOTpercent']

            rightChoice = np.append(rightChoice, rightChoicethisMouse)
            #valid = np.append(valid, validThisMouse)
            targetFrequency = np.append(targetFrequency, targetFrequencythisMouse)
            valid = np.ones(len(rightChoice), dtype=bool)

        (possibleValues, fractionHitsEachValue, ciHitsEachValue, nTrialsEachValue, nHitsEachValue) = behavioranalysis.calculate_psychometric(rightChoice, targetFrequency, valid)


        #curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues, fractionHitsEachValue)
        if thisTask == 'FT':
            possibleValuesFt = possibleValues
            fractionHitsEachValueFt[indMouse,:] = fractionHitsEachValue
            ciHitsEachValueFt[indMouse,:] = ciHitsEachValue
            nTrialsEachValueFt[indMouse,:] = nTrialsEachValue
            nHitsEachValueFt[indMouse,:] = nHitsEachValue

        elif thisTask == 'VOT':
            possibleValuesVot = possibleValues
            fractionHitsEachValueVot[indMouse,:] = fractionHitsEachValue
            ciHitsEachValueVot[indMouse,:] = ciHitsEachValue
            nTrialsEachValueVot[indMouse,:] = nTrialsEachValue
            nHitsEachValueVot[indMouse,:] = nHitsEachValue


fractionHitsEachValueFtAvg = np.mean(fractionHitsEachValueFt, axis=0)
ciHitsEachValueFtAvg = np.mean(ciHitsEachValueFt, axis=0)
nTrialsEachValueFtAvg = np.mean(nTrialsEachValueFt, axis=0)
nHitsEachValueFtAvg = np.mean(nHitsEachValueFt, axis=0)


fractionHitsEachValueVotAvg = np.mean(fractionHitsEachValueVot, axis=0)
ciHitsEachValueVotAvg = np.mean(ciHitsEachValueVot, axis=0)
nTrialsEachValueVotAvg = np.mean(nTrialsEachValueVot, axis=0)
nHitsEachValueVotAvg = np.mean(nHitsEachValueVot, axis=0)


curveParamsFt, pcovFt = scipy.optimize.curve_fit(extrastats.psychfun, possibleValuesFt, fractionHitsEachValueFtAvg)
curveParamsVot, pcovVot = scipy.optimize.curve_fit(extrastats.psychfun, possibleValuesVot, fractionHitsEachValueVotAvg)


if 1:

    xPadFt = 0.2 * (possibleValuesFt[-1] - possibleValuesFt[0])
    fitxvalFt = np.linspace(possibleValuesFt[0]-xPadFt, possibleValuesFt[-1]+xPadFt, 40)
    fityvalFt = extrastats.psychfun(fitxvalFt, *curveParamsFt)
    xTicks = np.arange(-1, 1.5, 0.5)
    fontSizeLabels = 12
    hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
    semBarsFT = np.array((np.std(fractionHitsEachValueFt,0)/np.sqrt(nFTmice)))
    plt.errorbar(possibleValuesFt, fractionHitsEachValueFtAvg*100, semBarsFT*100, marker = 'o', ms = 4, linestyle = '', color = 'k')
    #semBarsFT = np.array([(fractionHitsEachValueFtAvg - semBars), (semBars + fractionHitsEachValueFtAvg)])
    #(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesFt, fractionHitsEachValueFtAvg, semBars, xTicks=None, xscale='linear')
    #pline.set_visible(False)
    plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
    plt.xlabel('FT slope (A.U.)', fontsize=fontSizeLabels)
    plt.title('FT cohort average (n=3)', fontsize=fontSizeLabels, fontweight='bold')
    plt.grid(True, axis='y', color='0.9')

    for indMouse, thisMouse in enumerate(studyparams.BEHAVIOR_MICE['FT']):
        curveParamsFt_oneMouse, pcovFt_oneMouse = scipy.optimize.curve_fit(extrastats.psychfun, possibleValuesFt, fractionHitsEachValueFt[indMouse])
        xPadFt = 0.2 * (possibleValuesFt[-1] - possibleValuesFt[0])
        fitxvalFt = np.linspace(possibleValuesFt[0]-xPadFt, possibleValuesFt[-1]+xPadFt, 40)
        fityvalFt = extrastats.psychfun(fitxvalFt, *curveParamsFt_oneMouse)
        xTicks = np.arange(-1, 1.5, 0.5)
        fontSizeLabels = 12
        hfit = plt.plot(fitxvalFt, 100*fityvalFt, '--', linewidth=1, alpha = 0.6)
        #(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesFt, fractionHitsEachValueFt[indMouse], ciHitsEachValueFt[indMouse], xTicks=None, xscale='linear')
        #pline.set_visible(False)
        #plt.setp(pdots, ms = 4, mec = 'k', mfc = 'none')

    plt.show()

    plt.figure()
    xPadVot = 0.2 * (possibleValuesVot[-1] - possibleValuesVot[0])
    fitxvalVot = np.linspace(possibleValuesVot[0]-xPadVot, possibleValuesVot[-1]+xPadVot, 40)
    fityvalVot = extrastats.psychfun(fitxvalVot, *curveParamsVot)
    fontSizeLabels = 12
    hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
    semBarsVOT = np.array((np.std(fractionHitsEachValueVot,0)/np.sqrt(nVOTmice)))
    plt.errorbar(possibleValuesVot, fractionHitsEachValueVotAvg*100, semBarsVOT*100, marker = 'o', ms = 4, linestyle = '', color = 'k')

    #semBars = np.array((np.std(fractionHitsEachValueVot,0)/np.sqrt(nVOTmice)))
    #semBarsVOT = np.array([(fractionHitsEachValueVotAvg - semBars), (semBars + fractionHitsEachValueVotAvg)])
    #(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesVot, fractionHitsEachValueVotAvg, semBarsVOT, xTicks=None, xscale='linear')
    #pline.set_visible(False)
    plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
    plt.xlabel('Vot (A.U.)', fontsize=fontSizeLabels)
    plt.title('VOT cohort average (n=9)', fontsize=fontSizeLabels, fontweight='bold')
    plt.grid(True, axis='y', color='0.9')

    for indMouse, thisMouse in enumerate(studyparams.BEHAVIOR_MICE['VOT']):
        curveParamsVot_oneMouse, pcovVot_oneMouse = scipy.optimize.curve_fit(extrastats.psychfun, possibleValuesVot, fractionHitsEachValueVot[indMouse])

        fitxvalVot = np.linspace(possibleValuesVot[0]-xPadVot, possibleValuesVot[-1]+xPadVot, 40)
        fityvalVot = extrastats.psychfun(fitxvalVot, *curveParamsVot_oneMouse)
        xTicks = np.arange(-1, 1.5, 0.5)
        fontSizeLabels = 12
        hfit = plt.plot(fitxvalVot, 100*fityvalVot, '--', linewidth=1, alpha = 0.6)
        #plt.errorbar(possibleValuesVot, fractionHitsEachValueVot[indMouse]*100, yerr = ciHitsEachValueVot[indMouse]*100, linestyle = '', marker = 'o', ms = 4, alpha = 0.6)
        #(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValuesVot, fractionHitsEachValueVot[indMouse], ciHitsEachValueVot[indMouse], xTicks=None, xscale='linear')
        #pdots = plt.plot(possibleValuesVot, 100*fractionHitsEachValueVot[indMouse], 'o', ms=4, alpha = 0.6)
        #(pline, pcaps, pbars) = plt.errorbar(possibleValuesVot, fractionHitsEachValueVot[indMouse]*100, yerr = (ciHitsEachValueVot[indMouse] - fractionHitsEachValueVot[indMouse])*100, marker = 'o', mec = 'k', ms = 4, elinewidth = 2)
        #pline.set_visible(False)
        #plt.setp(pbars, lw = 2)
        #plt.setp(pdots, ms = 4, mec = 'k', mfc = 'none')
        #pdots = plt.plot(possibleValues, 100*fractionHitsEachValue, 'o', mec='none', mfc='k', ms=8)

    plt.show()

np.savez(figDataFullPath, possibleValuesFt = possibleValuesFt, fractionHitsEachValueFt = fractionHitsEachValueFt, fractionHitsEachValueFtAvg = fractionHitsEachValueFtAvg, semBarsFT = semBarsFT, ciHitsEachValueFt = ciHitsEachValueFt, nTrialsEachValueFtAvg = nTrialsEachValueFtAvg, nHitsEachValueFtAvg = nHitsEachValueFtAvg, curveParamsFt = curveParamsFt, pcovFt = pcovFt, possibleValuesVot = possibleValuesVot, fractionHitsEachValueVot = fractionHitsEachValueVot,  fractionHitsEachValueVotAvg = fractionHitsEachValueVotAvg, semBarsVOT = semBarsVOT, ciHitsEachValueVot = ciHitsEachValueVot, nTrialsEachValueVotAvg = nTrialsEachValueVotAvg, nHitsEachValueVotAvg = nHitsEachValueVotAvg, curveParamsVot = curveParamsVot, pcovVot = pcovVot)
print('saved to ' f'{figDataFullPath}')
