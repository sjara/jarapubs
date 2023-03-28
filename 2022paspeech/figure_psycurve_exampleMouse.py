'''
This figure plots 4 psychometric curves: 1 example mouse and cohort averages for both the FT and the VOT cohorts.
'''

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import scipy.stats as stats
import figparams
import studyparams
from importlib import reload
reload(figparams)

FIGNAME = 'figure_behavior'
figDataFileExample = 'data_example_mice_psycurves.npz'
figDataFileCohort = 'data_cohort_average_psycurves.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPathExample = os.path.join(figDataDir, figDataFileExample)
figDataFullPathCohort = os.path.join(figDataDir, figDataFileCohort)

PANELS = [2,2] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_psycurves' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [9, 7] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.54]   # Horiz position for panel labels
labelPosY = [0.95, 0.45]    # Vert position for panel labels

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
gsMain.update(left=0.15, right=0.98, top=0.92, bottom=0.08, wspace=0.25, hspace=0.4)

# -- Panel: name of panel --

# -- Panel: FT example mouse --
ax1 = plt.subplot(gsMain[0, 0])
ax1.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
ax1.annotate('n = 1 mouse \n 8 sessions', xy = (0,75), xycoords = 'data')
xPadFt = 0.2 * (figData['possibleValuesFt'][-1] - figData['possibleValuesFt'][0])
fitxvalFt = np.linspace(figData['possibleValuesFt'][0] - xPadFt, figData['possibleValuesFt'][-1] + xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *figData['curveParamsFt'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','da'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesFt'], figData['fractionHitsEachValueFt'], figData['ciHitsEachValueFt'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.xticks([0,20, 40, 60, 80, 100],['9', '6', '2', '-2', '-6', '-9'])
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Formant Transition slope (oct/s)', fontsize=fontSizeLabels)
plt.title('Example mouse: FT', fontsize=fontSizeLabels, fontweight='bold', pad=10)
plt.grid(True, axis='y', color='0.9')
## ADD ANNOTATE FOR /BA/ <-> /DA/ etc

if PANELS[0]:
    # Plot stuff
    pass

## -- Panel: VOT example mouse --
ax2 = plt.subplot(gsMain[1,0])
ax2.annotate('C', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
ax2.annotate('n = 1 mouse \n 8 sessions', xy = (0,75), xycoords = 'data')
#ax2.annotate('/ba/', xy = (0, 103), xycoords = 'data', annotation_clip = 0)
#ax2.annotate('/pa/', xy = (100, 103), xycoords = 'data', annotation_clip = 0)
#ax2.annotate(xy = (50, 103), xycoords = 'data', annotation_clip = 0, arrowprops = {arrowstyle: ``'<|-|>'``},)
xPadVot = 0.2 * (figData['possibleValuesVot'][-1] - figData['possibleValuesVot'][0])
fitxvalVot = np.linspace(figData['possibleValuesVot'][0] - xPadVot, figData['possibleValuesVot'][-1] + xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *figData['curveParamsVot'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','pa'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesVot'], figData['fractionHitsEachValueVot'], figData['ciHitsEachValueVot'], xTicks = None, xscale='linear')
pline.set_visible(False)
plt.xticks([0,20,40,60,80,100], ['2', '4', '8', '16', '32', '64'])
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Voice Onset Time (ms)', fontsize=fontSizeLabels)
plt.title('Example mouse: VOT', fontsize=fontSizeLabels, fontweight='bold', pad=10)
plt.grid(True, axis='y', color='0.9')


# -- Panel: FT cohort average --
ax3 = plt.subplot(gsMain[0, 1])
ax3.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
ax3.annotate('n = 3 mice \n 8 sessions', xy = (0,75), xycoords = 'data')

xPadFt = 0.2 * (cohortData['possibleValuesFt'][-1] - cohortData['possibleValuesFt'][0])
fitxvalFt = np.linspace(cohortData['possibleValuesFt'][0]-xPadFt, cohortData['possibleValuesFt'][-1]+xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *cohortData['curveParamsFt'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','da'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesFt'], cohortData['fractionHitsEachValueFtAvg'], cohortData['semBarsFT'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.xticks([0,20, 40, 60, 80, 100],['9', '6', '2', '-2', '-6', '-9'])
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Formant Transisiton slope (oct/s)', fontsize=fontSizeLabels)
plt.title('Cohort average: FT', fontsize=fontSizeLabels, fontweight='bold', pad=10)
plt.grid(True, axis='y', color='0.9')

## -- Panel: VOT cohort average --
ax4 = plt.subplot(gsMain[1,1])
ax4.annotate('D', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
ax4.annotate('n = 9 mice \n 8 sessions', xy = (0,75), xycoords = 'data')
xPadVot = 0.2 * (cohortData['possibleValuesVot'][-1] - cohortData['possibleValuesVot'][0])
fitxvalVot = np.linspace(cohortData['possibleValuesVot'][0]-xPadVot, cohortData['possibleValuesVot'][-1]+xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *cohortData['curveParamsVot'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','pa'])
fontSizeLabels = 12
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesVot'], cohortData['fractionHitsEachValueVotAvg'], cohortData['semBarsVOT'], xTicks=None, xscale='linear')
pline.set_visible(False)
plt.xticks([0,20,40,60,80,100], ['2', '4', '8', '16', '32', '64'])
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Voice Onset Time (ms)', fontsize=fontSizeLabels)
plt.title('Cohort average: VOT', fontsize=fontSizeLabels, fontweight='bold', pad=10)
plt.grid(True, axis='y', color='0.9')

if PANELS[1]:
    # Plot stuff
    pass

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
