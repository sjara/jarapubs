"""
For figure showing example performance for one mouse.
"""

import os
import sys
import numpy as np
import scipy
import studyparams
import studyutils
from jaratoolbox import settings
from jaratoolbox import loadbehavior
from jaratoolbox import behavioranalysis
from jaratoolbox import extrastats

SAVE_DATA = 0

FIGNAME = 'example_performance'
figDataFile = 'example_performance_{}.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

# -- Optional --
if not os.path.exists(figDataDir):
    #os.mkdir(figDataDir)
    print('Please create folder: {}'.format(figDataDir)); sys.exit()

subject = 'pamo009'

# -- Load learning data (stage 3) --
dframe = studyutils.load_stage3(ndays=26, excludeAntibias=True)
perfEachDayStage3 = dframe.loc[subject]

# -- Load psychometric data (stage 4) --
#session = '20220305a' # pamo010 Nicest
#session = '20220308a' # pamo009
session = '20220401a' # pamo009

paradigm = '2afc'
varlist = ['outcome', 'valid', 'choice', 'targetFMslope', 'startFreq', 'endFreq', 'targetDuration']
behavFile = loadbehavior.path_to_behavior_data(subject, paradigm, session)
bdata = loadbehavior.BehaviorData(behavFile, varlist)

valid = bdata['valid'].astype(bool)
leftwardChoice = bdata['choice']==bdata.labels['choice']['left']

#targetParamValue = bdata['targetFMslope']
#possibleParamValue = np.unique(targetParamValue)
#nParamValues = len(possibleParamValue)
# -- The FM slope is in units of octaves/sec --
targetDuration = bdata['targetDuration'][-1]
targetParamValue = np.log2(bdata['endFreq']/bdata['startFreq'])/targetDuration

(possibleValues, fractionLeftEachValue, ciLeftEachValue, nTrialsEachValue, nLeftEachValue)=\
    behavioranalysis.calculate_psychometric(leftwardChoice, targetParamValue, valid)

# -- Fit sigmoidal curve --
par0 = [0, 0.5, 0, 0]
bounds = [[-np.inf, 0.08, 0, 0], [np.inf, np.inf, 0.5, 0.5]]
curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleValues,
                                              fractionLeftEachValue, p0=par0, bounds=bounds)

if 1:
    import matplotlib.pyplot as plt
    from jaratoolbox import extraplots
    plt.clf()
    xPad = 0.2 * (possibleValues[-1] - possibleValues[0])
    fitxval = np.linspace(possibleValues[0]-xPad, possibleValues[-1]+xPad, 40)
    fityval = extrastats.psychfun(fitxval, *curveParams)
    xTicks = np.arange(-6, 7, 2)
    hfit = plt.plot(fitxval, 100*fityval, '-', lw=3, color='k')
    (pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValues,
                                                                fractionLeftEachValue,
                                                                ciLeftEachValue, xTicks=xTicks,
                                                                xscale='linear')
    pline.set_visible(False)
    plt.xlim([-6, 6])
    plt.ylabel('Leftward choice (%)')
    plt.xlabel('FM slope (oct/s)')
    extraplots.boxoff(plt.gca())
    plt.show()
    

if SAVE_DATA:
    outputFile = figDataFullPath.format(subject)
    np.savez(outputFile, script=scriptFullPath, subject=subject, perfEachDayStage3=perfEachDayStage3,
             possibleFMslopes=possibleValues, fractionLeftEachValue=fractionLeftEachValue,
             ciLeftEachValue=ciLeftEachValue, curveParams=curveParams, psycurveSession=session)
    print(f'Saved {outputFile}')

