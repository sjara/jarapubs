'''
This script provides a template for scripts to generate figures.
Replace this comment with a description of which figure will be plotted by this script.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
#import scipy.stats as stats
import figparams
import studyparams

FIGNAME = 'figure_behavior'
figDataFileExample = 'data_example_mice_psycurves.npz'
figDataFileCohort = 'data_cohort_average_psycurves.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPathExample = os.path.join(figDataDir, figDataFileExample)
figDataFullPathCohort = os.path.join(figDataDir, figDataFileCohort)

PANELS = [2,2] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_behavior' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.9, 0.48]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
sound = figparams.colors['sound']
ftCond = figparams.colors['ft']
votCond = figparams.colors['vot']


# -- Load data --
figData = np.load(figDataFullPathExample)
cohortData = np.load(figDataFullPathCohort)

# -- Processed data --
pass

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 2)
gsMain.update(left=0.15, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)

# -- Panel: name of panel --

# -- Panel: FT example mouse --
ax1 = plt.subplot(gsMain[0, 0])
#ax1.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
xPadFt = 0.2 * (figData['possibleValuesFt'][-1] - figData['possibleValuesFt'][0])
fitxvalFt = np.linspace(figData['possibleValuesFt'][0] - xPadFt, figData['possibleValuesFt'][-1] + xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *figData['curveParamsFt'])
xTicks = np.arange(-1, 1.5, 0.5)
fontSizeLabels = 12
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesFt'], figData['fractionHitsEachValueFt'], figData['ciHitsEachValueFt'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('FT slope (A.U.)', fontsize=fontSizeLabels)
plt.title('FT example mouse', fontsize=fontSizeLabels, fontweight='bold')
plt.grid(True, axis='y', color='0.9')


if PANELS[0]:
    # Plot stuff
    pass

## -- Panel: VOT example mouse --
ax2 = plt.subplot(gsMain[0,1])
#ax2.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
xPadVot = 0.2 * (figData['possibleValuesVot'][-1] - figData['possibleValuesVot'][0])
fitxvalVot = np.linspace(figData['possibleValuesVot'][0] - xPadVot, figData['possibleValuesVot'][-1] + xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *figData['curveParamsVot'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesVot'], figData['fractionHitsEachValueVot'], figData['ciHitsEachValueVot'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Vot (A.U.)', fontsize=fontSizeLabels)
plt.title('VOT example mouse', fontsize=fontSizeLabels, fontweight='bold')
plt.grid(True, axis='y', color='0.9')


# -- Panel: FT cohort average --
ax3 = plt.subplot(gsMain[1, 0])
xPadFt = 0.2 * (cohortData['possibleValuesFt'][-1] - cohortData['possibleValuesFt'][0])
fitxvalFt = np.linspace(cohortData['possibleValuesFt'][0]-xPadFt, cohortData['possibleValuesFt'][-1]+xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *cohortData['curveParamsFt'])
xTicks = np.arange(-1, 1.5, 0.5)
fontSizeLabels = 12
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesFt'], cohortData['fractionHitsEachValueFtAvg'], cohortData['ciHitsEachValueFtAvg'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('FT slope (A.U.)', fontsize=fontSizeLabels)
plt.title('FT cohort average (n=3)', fontsize=fontSizeLabels, fontweight='bold')
plt.grid(True, axis='y', color='0.9')

## -- Panel: VOT cohort average --
ax4 = plt.subplot(gsMain[1,1])
xPadVot = 0.2 * (cohortData['possibleValuesVot'][-1] - cohortData['possibleValuesVot'][0])
fitxvalVot = np.linspace(cohortData['possibleValuesVot'][0]-xPadVot, cohortData['possibleValuesVot'][-1]+xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *cohortData['curveParamsVot'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesVot'], cohortData['fractionHitsEachValueVotAvg'], cohortData['ciHitsEachValueVotAvg'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Vot (A.U.)', fontsize=fontSizeLabels)
plt.title('VOT cohort average (n=9)', fontsize=fontSizeLabels, fontweight='bold')
plt.grid(True, axis='y', color='0.9')

if PANELS[1]:
    # Plot stuff
    pass

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
