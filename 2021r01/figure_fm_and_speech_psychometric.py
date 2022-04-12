"""
Based on ~/src/jaratest/common/2020nsf/figure_example_speech_discrimination.py

"""


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from jaratoolbox import behavioranalysis
#from jaratoolbox import loadbehavior
from jaratoolbox import extraplots
from statsmodels.stats.proportion import proportion_confint #Used to compute confidence interval for the error bars. 
from jaratoolbox import settings 
#import figparams

SAVE_FIGURE = 1
outputDir = '/tmp/'

figFilename = 'fm_and_speech_psycurves'
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [2.6, 5]

fontSizeLabels = 16
fontSizeTicks = 16

plotColor = 'k'

# Spectral bili007 from 2018-10-24 through 2018-10-28
# Temporal bili006 from 2018-10-26 through 2018-11-02

    
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')
gs = gridspec.GridSpec(2, 1)
gs.update(top=0.96, left=0.4, right=0.95, bottom=0.12, wspace=0.55, hspace=0.5)

# ---------------- FM ----------------
subject = 'frem002'
#sessions = ['20211012a','20211013a','20211014a','20211015a',
#            '20211016a','20211017a','20211018a','20211019a']
sessions = ['20211012a','20211016a','20211017a','20211018a']
#sessions = ['20211019a'] # 12 16 17 18
paradigm = '2afc'

bdata = behavioranalysis.load_many_sessions(subject,paradigm=paradigm,sessions=sessions)
validTrials = bdata['valid'].astype(bool)
correctTrials = bdata['outcome']==bdata.labels['outcome']['correct']
rightwardChoice = bdata['choice']==bdata.labels['choice']['right']

targetParamValue = bdata['targetFMslope']
possibleParamValue = np.unique(targetParamValue)
nParamValues = len(possibleParamValue)

# -- Calculate max FM slope values --
highFreq = bdata['highFreq'][-1]
lowFreq = bdata['lowFreq'][-1]
targetDuration = bdata['targetDuration'][-1]
maxFMslope = np.log2(highFreq/lowFreq)/targetDuration
print(f'Max FM slope: {maxFMslope:0.2f} oct/s')

(possibleValues, fractionHitsEachValue, ciHitsEachValue, nTrialsEachValue, nHitsEachValue)=\
    behavioranalysis.calculate_psychometric(rightwardChoice, targetParamValue, validTrials)

ax0 = plt.subplot(gs[0, 0])
xTicks = np.arange(-1, 1.2, 0.2)

# -- Flip for display purposes --
fractionHitsEachValue = np.flip(fractionHitsEachValue)
ciHitsEachValue = np.fliplr(ciHitsEachValue)

(pline, pcaps, pbars, pdots) = \
    extraplots.plot_psychometric(possibleValues, fractionHitsEachValue,
                                 ciHitsEachValue, xTicks=xTicks, xscale='linear')
plt.setp([pline, pcaps, pbars, pdots], color=plotColor)
plt.setp(pdots, mfc=plotColor, mec='none')
#plt.axhline(y=50, color='0.85', ls='--')
#plt.minorticks_off()
plt.ylim([-5,105])
#plt.ylabel('Rightward trials (%)', fontsize=fontSizeLabels)
plt.ylabel('Reported\ndown (%)', fontsize=fontSizeLabels)
plt.xlabel('FM slope (oct/s)', fontsize=fontSizeLabels)
sessionStr = '{} - {}'.format(sessions[0],sessions[-1])
titleStr = f'{subject}: {sessionStr}'
#plt.title(titleStr, fontsize=fontSizeLabels, fontweight='bold')
ax0.set_yticks([0, 50, 100])
ax0.set_xticks([-1, 0, 1])
ax0.set_xticklabels(['6','0','-6'])
extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
extraplots.boxoff(ax0)

# ---------------- SPEECH ----------------
CASE = 0
if CASE==0:
    subject = 'bili007'
    sessions = ['20181024a','20181025a','20181027a','20181028a']
    taskLabel = 'spectral'
    #xLabel = 'Formant transition slope (Hz/oct)'
    #xLabel = 'Formant transition slope (%)'
    xLabel = 'FT slope (%)'
elif CASE==1:
    subject = 'bili006'
    sessions = ['20181027a','20181029a','20181030a','20181101a']
    taskLabel = 'temporal'
    xLabel = 'Voice Onset Time (ms)'

bdata = behavioranalysis.load_many_sessions(subject,sessions)

targetPercentage = bdata['targetFrequency'] # I used name 'frequency' initially
choiceRight = bdata['choice']==bdata.labels['choice']['right']
valid=bdata['valid']& (bdata['choice']!=bdata.labels['choice']['none'])
(possibleValues,fractionHitsEachValue,ciHitsEachValue,nTrialsEachValue,nHitsEachValue)=\
       behavioranalysis.calculate_psychometric(choiceRight,targetPercentage,valid)

ax1 = plt.subplot(gs[1, 0])
#plt.title('{0} [{1}-{2}]'.format(subject,sessions[0],sessions[-1]))

(pline, pcaps, pbars, pdots) = extraplots.plot_psychometric(possibleValues,fractionHitsEachValue,
                                                            ciHitsEachValue, xscale='linear')

plt.setp([pline,pcaps,pbars,pdots], color=plotColor)
plt.setp(pdots, mfc=plotColor, mec='none')
ax1.set_yticks([0, 50, 100])
ax1.set_xticks([0, 50, 100])
ax1.set_xticklabels(['-100','','100'])
plt.xlabel(xLabel,fontsize=fontSizeLabels)
#plt.ylabel('Rightward trials (%)',fontsize=fontSizeLabels)
plt.ylabel('Reported\n/da/ (%)',fontsize=fontSizeLabels)
extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
extraplots.boxoff(ax1)



plt.show()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
