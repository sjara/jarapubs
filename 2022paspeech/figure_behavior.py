'''
This creates figure 1 for 2022paspeech:
 A. plots 3 spectrograms: /ba/, /pa/, and /da/ at 8x human freq range
 B. 2AFC cartoon
 C - F. Psycurves: 1 example mouse and cohort averages for both the FT and the VOT cohorts.
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
from jaratoolbox import soundanalysis
from jaratoolbox import colorpalette as cp
import scipy.stats as stats
from scipy.io import wavfile
import scipy.optimize
import figparams
import studyparams
from matplotlib.font_manager import findfont, FontProperties
from importlib import reload
font = findfont(FontProperties(family = ['Helvetica']))
reload(figparams)

FIGNAME = 'figure_behavior'
SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_behavior' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 5.5] # In inches

PANELS = [3,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
colorAllFits = cp.TangoPalette['Aluminium3']

labelPosX = [0.04, 0.42, 0.72] # Horiz position for panel labels
labelPosY = [0.95, 0.67, 0.47]    # Vert position for panel labels

plt.figure()
#gsMain = gridspec.GridSpec(2, 1, height_ratios=[0.3, 0.7])
gsMain = gridspec.GridSpec(1, 2, width_ratios = [0.35, 0.65])
gsMain.update(left=0.1, right=0.95, top=0.89, bottom=0.12, wspace=0.45, hspace=0.4)

gsLeft = gsMain[0].subgridspec(3, 2, height_ratios = [0.3, 0.3, 0.3], hspace = 0.6, wspace = 0.3)
axTask = plt.subplot(gsLeft[0,0])
axDa = plt.subplot(gsLeft[2, 0])
axBa = plt.subplot(gsLeft[1, 0])
axPa = plt.subplot(gsLeft[1, 1])
axCbar = plt.subplot(gsLeft[2,1])
plt.subplots_adjust(top = 0.95, bottom = 0.1, hspace = 0.4, left = 0.05)


gsPsy = gsMain[1].subgridspec(2,2, wspace = 0.3, hspace = 0.6)
axExFt = plt.subplot(gsPsy[1,0])
axExVot = plt.subplot(gsPsy[0,0])
axPopFt = plt.subplot(gsPsy[1,1])
axPopVot = plt.subplot(gsPsy[0,1])
plt.subplots_adjust(top = 0.9, bottom = 0.1, hspace = 0.4, left = 0.05)


axTask.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axBa.annotate('B', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axExVot.annotate('C', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axPopVot.annotate('D', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axExFt.annotate('E', xy=(labelPosX[1],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axPopFt.annotate('F', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')





# -- Panel A: spectrograms /ba/, /da/, /pa/ --
soundsDir = 'H:\\jarasounds\\ft_vot_8x_20220115' #HARDCODED

# -- /da/ --
#ax1 = plt.subplot(gsMain[0,0])
thisSound = f'syllable_8x_vot000_ft100.wav'

plt.sca(axDa)
soundFile = os.path.join(soundsDir, thisSound)
samplingRate, wave = wavfile.read(soundFile)

[sgramF, sgramT, sgramV] = scipy.signal.stft(wave, fs=samplingRate,
                               window='hanning', nperseg=2048, noverlap=1024)
INTERP = 'nearest'
sgramVsq = np.abs(sgramV**2)
minVal =  np.min(sgramVsq[sgramVsq!=0])
sgramVsq[sgramVsq==0] = minVal
sgramVals = np.log10(sgramVsq)

intensityRange = sgramVals.max()-sgramVals.min()
VMAX=None; VMIN = sgramVals.min()+0.25*intensityRange
#plt.clf()
plt.imshow(sgramVals, cmap='viridis', aspect='auto',
           interpolation=INTERP, vmin=VMIN, vmax=VMAX,
           extent=(sgramT[0],sgramT[-1],sgramF[-1],sgramF[0]))
plt.gca().invert_yaxis()
plt.ylabel('Frequency (kHz)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.yticks(np.arange(0, 28000, 5000), ['0', '5', '10', '15', '20', '25'], fontsize=fontSizeTicks)
plt.ylim([0,28000])
plt.xlim([0,.1])
plt.xticks(np.arange(0, .11, 0.05),['0', '50', '100'], fontsize=fontSizeTicks)
plt.xlabel('Time (ms)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.title('/da/', fontsize=fontSizeLabels, fontweight='bold')

cbar = plt.colorbar(ax = axCbar, orientation = 'vertical', location = 'left', ticks = [], pad = 0.4) #shrink = 0.8,
axCbar.annotate('Intensity', xy = (0.26, 0.16), xycoords='subfigure fraction', fontsize= fontSizeLabels, rotation = 270)
axCbar.set_axis_off()

# -- /ba/ --
#ax2 = plt.subplot(gsMain[0,1])
plt.sca(axBa)
thisSound = f'syllable_8x_vot000_ft000.wav'
soundFile = os.path.join(soundsDir, thisSound)
samplingRate, wave = wavfile.read(soundFile)
#soundanalysis.plot_spectrogram(wave, samplingRate)

[sgramF, sgramT, sgramV] = scipy.signal.stft(wave, fs=samplingRate,
                               window='hanning', nperseg=2048, noverlap=1024)
INTERP = 'nearest'
sgramVsq = np.abs(sgramV**2)
minVal =  np.min(sgramVsq[sgramVsq!=0])
sgramVsq[sgramVsq==0] = minVal
sgramVals = np.log10(sgramVsq)

intensityRange = sgramVals.max()-sgramVals.min()
VMAX=None; VMIN = sgramVals.min()+0.25*intensityRange
#plt.clf()
plt.imshow(sgramVals, cmap='viridis', aspect='auto',
           interpolation=INTERP, vmin=VMIN, vmax=VMAX,
           extent=(sgramT[0],sgramT[-1],sgramF[-1],sgramF[0]))
#plt.ylim(np.array([200,0]))
plt.gca().invert_yaxis()
plt.ylabel('Frequency (kHz)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.yticks(np.arange(0, 28000, 5000), ['0', '5', '10', '15', '20', '25'], fontsize=fontSizeTicks)
plt.ylim([0,28000]) ##CHANGE THIS TO SHARE AX W/DA
plt.xlim([0,.1])
plt.xticks(np.arange(0, .11, 0.05), labels = '', fontsize=fontSizeTicks)
#plt.xlabel('Time (ms)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.title('/ba/', fontsize=fontSizeLabels, fontweight='bold')
plt.xlabel('')

# -- /pa/ --
#ax3 = plt.subplot(gsMain[0,2])
plt.sca(axPa)
thisSound = f'syllable_8x_vot100_ft000.wav'
soundFile = os.path.join(soundsDir, thisSound)
samplingRate, wave = wavfile.read(soundFile)
#soundanalysis.plot_spectrogram(wave, samplingRate)

[sgramF, sgramT, sgramV] = scipy.signal.stft(wave, fs=samplingRate,
                               window='hanning', nperseg=2048, noverlap=1024)
INTERP = 'nearest'
sgramVsq = np.abs(sgramV**2)
minVal =  np.min(sgramVsq[sgramVsq!=0])
sgramVsq[sgramVsq==0] = minVal
sgramVals = np.log10(sgramVsq)

intensityRange = sgramVals.max()-sgramVals.min()
VMAX=None; VMIN = sgramVals.min()+0.25*intensityRange
#plt.clf()
plt.imshow(sgramVals, cmap='viridis', aspect='auto',
           interpolation=INTERP, vmin=VMIN, vmax=VMAX,
           extent=(sgramT[0],sgramT[-1],sgramF[-1],sgramF[0]))
#plt.ylim(np.array([200,0]))
plt.gca().invert_yaxis()
plt.ylim([0,28000])
plt.xlim([0,.1])
plt.xticks(np.arange(0, .11, 0.05),['0', '50', '100'], fontsize=fontSizeTicks)
plt.xlabel('Time (ms)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.title('/pa/', fontsize=fontSizeLabels, fontweight='bold')
plt.yticks(np.arange(0, 28000, 5000), labels = '', fontsize=fontSizeTicks)
plt.ylabel('')


# -- Panel B: 2AFC task cartoon --
#ax4 = plt.subplot(gsMain[0,3:])
plt.sca(axTask)
axTask.set_axis_off()

# -- Panel C - F: Psycurves --
figDataFileExample = 'data_example_mice_psycurves.npz'
figDataFileCohort = 'data_cohort_average_psycurves.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPathExample = os.path.join(figDataDir, figDataFileExample)
figDataFullPathCohort = os.path.join(figDataDir, figDataFileCohort)


labelPosX = [0.07, 0.54, 0.7]   # Horiz position for panel labels
labelPosY = [0.95, 0.45]    # Vert position for panel labels

# -- Load data --
figData = np.load(figDataFullPathExample)
cohortData = np.load(figDataFullPathCohort)

# -- Plot results --
#fig = plt.gcf()
#fig.clf()
#fig.set_facecolor('w')


# -- Panel C: FT example mouse --
#ax5 = plt.subplot(gsMain[1, 0:3])
plt.sca(axExFt)
#ax5.annotate('C', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axExFt.annotate('n = 1 mouse \n 8 sessions', xy = (0,75), xycoords = 'data', fontsize=fontSizeLabels)
xPadFt = 0.2 * (figData['possibleValuesFt'][-1] - figData['possibleValuesFt'][0])
fitxvalFt = np.linspace(figData['possibleValuesFt'][0] - xPadFt, figData['possibleValuesFt'][-1] + xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *figData['curveParamsFt'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','da'])
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesFt'], figData['fractionHitsEachValueFt'], figData['ciHitsEachValueFt'], xTicks=None, xscale='linear')
plt.setp(pdots, ms = 4)
pline.set_visible(False)
plt.xticks([0,20, 40, 60, 80, 100],['9', '6', '2', '-2', '-6', '-9'], fontsize=fontSizeTicks)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.xlabel('Formant Transition slope \n (oct/s)', fontsize=fontSizeLabels, labelpad = 0.3)
#plt.title('Example mouse: FT', fontsize=fontSizeLabels, pad=5)
#plt.grid(True, axis='y', color='0.9')
axExFt.spines["right"].set_visible(False)
axExFt.spines["top"].set_visible(False)

## -- Panel: VOT example mouse --
#ax6 = plt.subplot(gsMain[2,0:3])
plt.sca(axExVot)
#ax6.annotate('E', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axExVot.annotate('n = 1 mouse \n 8 sessions', xy = (0,75), xycoords = 'data')
#ax2.annotate('/ba/', xy = (0, 103), xycoords = 'data', annotation_clip = 0)
#ax2.annotate('/pa/', xy = (100, 103), xycoords = 'data', annotation_clip = 0)
#ax2.annotate(xy = (50, 103), xycoords = 'data', annotation_clip = 0, arrowprops = {arrowstyle: ``'<|-|>'``},)
xPadVot = 0.2 * (figData['possibleValuesVot'][-1] - figData['possibleValuesVot'][0])
fitxvalVot = np.linspace(figData['possibleValuesVot'][0] - xPadVot, figData['possibleValuesVot'][-1] + xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *figData['curveParamsVot'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','pa'])
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(figData['possibleValuesVot'], figData['fractionHitsEachValueVot'], figData['ciHitsEachValueVot'], xTicks = None, xscale='linear')
plt.setp(pdots, ms = 4)
pline.set_visible(False)
plt.xticks([0,20,40,60,80,100], ['2', '4', '8', '16', '32', '64'], fontsize=fontSizeTicks)
plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels, labelpad = 0.3)
plt.xlabel('Voice Onset Time (ms)', fontsize=fontSizeLabels, labelpad = 0.3)
#plt.title('Example mouse: VOT', fontsize=fontSizeLabels, pad=5)
#plt.grid(True, axis='y', color='0.9')
axExVot.spines["right"].set_visible(False)
axExVot.spines["top"].set_visible(False)


# -- Panel: FT cohort average --
#ax7 = plt.subplot(gsMain[1, 3:])
plt.sca(axPopFt)
#ax7.annotate('D', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axPopFt.annotate('n = 3 mice \n 8 sessions', xy = (0,75), xycoords = 'data')

xPadFt = 0.2 * (cohortData['possibleValuesFt'][-1] - cohortData['possibleValuesFt'][0])
fitxvalFt = np.linspace(cohortData['possibleValuesFt'][0]-xPadFt, cohortData['possibleValuesFt'][-1]+xPadFt, 40)
fityvalFt = extrastats.psychfun(fitxvalFt, *cohortData['curveParamsFt'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','da'])
hfit = plt.plot(fitxvalFt, 100*fityvalFt, '-', linewidth=2, color='k')
plt.errorbar(cohortData['possibleValuesFt'], cohortData['fractionHitsEachValueFtAvg']*100, cohortData['semBarsFT']*100, marker = 'o', linestyle = '', color = 'k', ms = 3)
#(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesFt'], cohortData['fractionHitsEachValueFtAvg'], cohortData['semBarsFT'], xTicks=None, xscale='linear')
#plt.setp(pdots, ms = 4)
#pline.set_visible(False)
for indMouse, thisMouse in enumerate(studyparams.BEHAVIOR_MICE['FT']):
    curveParamsFt_oneMouse, pcovFt_oneMouse = scipy.optimize.curve_fit(extrastats.psychfun, cohortData['possibleValuesFt'], cohortData['fractionHitsEachValueFt'][indMouse])
    xPadFt = 0.2 * (cohortData['possibleValuesFt'][-1] - cohortData['possibleValuesFt'][0])
    fitxvalFt = np.linspace(cohortData['possibleValuesFt'][0]-xPadFt, cohortData['possibleValuesFt'][-1]+xPadFt, 40)
    fityvalFt = extrastats.psychfun(fitxvalFt, *curveParamsFt_oneMouse)
    xTicks = np.arange(-1, 1.5, 0.5)
    hfit = plt.plot(fitxvalFt, 100*fityvalFt, c = colorAllFits, linewidth=1, alpha = 0.6)

valRange = cohortData['possibleValuesFt'][-1]-cohortData['possibleValuesFt'][0]
plt.xlim([cohortData['possibleValuesFt'][0]-0.1*valRange, cohortData['possibleValuesFt'][-1]+0.1*valRange])
plt.ylim([0,100])
plt.xticks([0,20, 40, 60, 80, 100],['9', '6', '2', '-2', '-6', '-9'], fontsize=fontSizeTicks)
#plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Formant Transisiton slope \n (oct/s)', fontsize=fontSizeLabels, labelpad = 0.3)
#plt.title('Cohort average: FT', fontsize=fontSizeLabels, pad=5)
#plt.grid(True, axis='y', color='0.9')
axPopFt.spines["right"].set_visible(False)
axPopFt.spines["top"].set_visible(False)

## -- Panel: VOT cohort average --
#ax8 = plt.subplot(gsMain[2,3:])
plt.sca(axPopVot)
#ax8.annotate('F', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

for indMouse, thisMouse in enumerate(studyparams.BEHAVIOR_MICE['VOT']):
    curveParamsVot_oneMouse, pcovVot_oneMouse = scipy.optimize.curve_fit(extrastats.psychfun, cohortData['possibleValuesVot'], cohortData['fractionHitsEachValueVot'][indMouse])
    xPadVot = 0.2 * (cohortData['possibleValuesVot'][-1] - cohortData['possibleValuesVot'][0])
    fitxvalVot = np.linspace(cohortData['possibleValuesVot'][0]-xPadVot, cohortData['possibleValuesVot'][-1]+xPadVot, 40)
    fityvalVot = extrastats.psychfun(fitxvalVot, *curveParamsVot_oneMouse)
    xTicks = np.arange(-1, 1.5, 0.5)
    hfit = plt.plot(fitxvalVot, 100*fityvalVot, c = colorAllFits, linewidth=1, alpha = 0.6)

valRange = cohortData['possibleValuesVot'][-1]-cohortData['possibleValuesVot'][0]
plt.xlim([cohortData['possibleValuesVot'][0]-0.1*valRange, cohortData['possibleValuesVot'][-1]+0.1*valRange])
plt.ylim([0,100])
plt.xticks([0,20,40,60,80,100], ['2', '4', '8', '16', '32', '64'], fontsize=fontSizeTicks)
#plt.ylabel('Rightward choice (%)', fontsize=fontSizeLabels)
plt.xlabel('Voice Onset Time (ms)', fontsize=fontSizeLabels, labelpad = 0.3)
#plt.title('Cohort average: VOT', fontsize=fontSizeLabels, pad=5)
#plt.grid(True, axis='y', color='0.9')

axPopVot.annotate('n = 9 mice \n 8 sessions', xy = (0,75), xycoords = 'data')
xPadVot = 0.2 * (cohortData['possibleValuesVot'][-1] - cohortData['possibleValuesVot'][0])
fitxvalVot = np.linspace(cohortData['possibleValuesVot'][0]-xPadVot, cohortData['possibleValuesVot'][-1]+xPadVot, 40)
fityvalVot = extrastats.psychfun(fitxvalVot, *cohortData['curveParamsVot'])
xTicks = (np.arange(-1, 1.5, 0.5),['ba','pa'])
hfit = plt.plot(fitxvalVot, 100*fityvalVot, '-', linewidth=2, color='k')
plt.errorbar(cohortData['possibleValuesVot'], cohortData['fractionHitsEachValueVotAvg']*100, cohortData['semBarsVOT']*100, marker = 'o', linestyle = '', color = 'k', ms = 3)
axPopVot.spines["right"].set_visible(False)
axPopVot.spines["top"].set_visible(False)
#(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(cohortData['possibleValuesVot'], cohortData['fractionHitsEachValueVotAvg'], cohortData['semBarsVOT'], xTicks=None, xscale='linear')
#plt.setp(pdots, ms = 4)
#pline.set_visible(False)

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
