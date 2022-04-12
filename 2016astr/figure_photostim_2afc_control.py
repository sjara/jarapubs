import os
import numpy as np
from matplotlib import pyplot as plt
from jaratoolbox import colorpalette as cp
from jaratoolbox import extraplots
from jaratoolbox import settings
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.lines as mlines
import scipy.stats as stats
import figparams
reload(figparams)

STUDY_NAME = '2016astr'
FIGNAME = 'photostim_2afc'
dataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)

PANELS = [1,1] # Which panels to plot

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'supp_photostim_2afc_control' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7,4]

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
labelDis = 0.1
labelPosX = [0.02, 0.5]   # Horiz position for panel labels
labelPosY = [0.95, 0.95]    # Vert position for panel labels

PHOTOSTIMCOLORS = {'no_laser':'k',
                   'laser_left':figparams.colp['stimLeft'],
                   'laser_right':figparams.colp['stimRight']}

SHAPES = {'d1pi014':'s',
          'd1pi015':'o',
          'd1pi016':'^'}

fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')
plt.hold(True)
gs = gridspec.GridSpec(1, 2)
gs.update(left=0.13, right=0.97, top=0.95, bottom=0.1, wspace=0.29, hspace=0.15)


# -- Panel A: L vs R hemi stim bias in d1pi and wildtype controls -- #
ax1 = plt.subplot(gs[0, 0])
ax1.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
# -- Load d1pi data -- #
summaryFilename = 'summary_photostim_percent_right_choice_change.npz'
summaryFullPath = os.path.join(dataDir,summaryFilename)
summary = np.load(summaryFullPath)

left014 = summary['d1pi014leftHemiStim']
left015 = summary['d1pi015leftHemiStim']
left016 = summary['d1pi016leftHemiStim']
right014 = summary['d1pi014rightHemiStim']
right015 = summary['d1pi015rightHemiStim']
right016 = summary['d1pi016rightHemiStim']


# -- To match figure 2: Select the first 10 sessions from each hemi each mouse so that they have equal number of sessions -- #
maxSessionNum = 10
left014 = left014[:maxSessionNum]
left015 = left015[:maxSessionNum]
left016 = left016[:maxSessionNum]
right014 = right014[:maxSessionNum]
right015 = right015[:maxSessionNum]
right016 = right016[:maxSessionNum]
print 'Only plotting the first {} sessions in the summary panel.'.format(maxSessionNum) 


# -- Load control data -- #
summaryFilename = 'summary_photostim_control_percent_right_choice_change.npz'
summaryFullPath = os.path.join(dataDir,summaryFilename)
summary = np.load(summaryFullPath)

left048 = summary['adap048leftHemiStim']
left056 = summary['adap056leftHemiStim']
right048 = summary['adap048rightHemiStim']
right056 = summary['adap056rightHemiStim']

ax1.axhline(y=0, color='k', linestyle='-')
np.random.seed(7) #2

for animal,leftData in zip(['d1pi014','d1pi015','d1pi016'],[left014,left015,left016]):
#for animal,leftData in zip(['d1pi014','d1pi015'],[left014,left015]):
    randOffset = 0.3*(np.random.rand(len(leftData))-0.5)
    ax1.plot(1+randOffset, 100*leftData, 'o', mec=PHOTOSTIMCOLORS['laser_left'], mfc='None')

for animal,controlLeftData in zip(['adap048','adap056'],[left048,left056]):
    randOffset = 0.3*(np.random.rand(len(controlLeftData))-0.5)
    ax1.plot(2+randOffset, 100*controlLeftData, 'o', mec=PHOTOSTIMCOLORS['laser_left'], mfc='None')

for animal,rightData in zip(['d1pi014','d1pi015','d1pi016'],[right014,right015,right016]):
#for animal,rightData in zip(['d1pi015','d1pi016'],[right015,right016]):
    randOffset = 0.3*(np.random.rand(len(rightData))-0.5)
    ax1.plot(3+randOffset, 100*rightData, 'o', mec=PHOTOSTIMCOLORS['laser_right'], mfc='None')

for animal,controlRightData in zip(['adap048','adap056'],[right048,right056]):
    randOffset = 0.3*(np.random.rand(len(controlRightData))-0.5)
    ax1.plot(4+randOffset, 100*controlRightData, 'o', mec=PHOTOSTIMCOLORS['laser_right'], mfc='None')


# -- Stats for summary panel in figure grouping all animals together -- #
leftStimChange = np.concatenate((left014,left015,left016))
rightStimChange = np.concatenate((right014,right015,right016))
#leftStimChange = np.concatenate((left014,left015))
#rightStimChange = np.concatenate((right015,right016))
meanLeftStim = np.mean(leftStimChange)
meanRightStim = np.mean(rightStimChange)
ax1.plot(0.3*np.array([-1,1])+1, 100*np.tile(meanLeftStim,2), lw=3, color=PHOTOSTIMCOLORS['laser_left'])
ax1.plot(0.3*np.array([-1,1])+3, 100*np.tile(meanRightStim,2), lw=3, color=PHOTOSTIMCOLORS['laser_right'])

controlLeftStimChange = np.concatenate((left048,left056))
controlRightStimChange = np.concatenate((right048,right056))
meanLeftStimCtrl = np.mean(controlLeftStimChange)
meanRightStimCtrl = np.mean(controlRightStimChange)
ax1.plot(0.3*np.array([-1,1])+2, 100*np.tile(meanLeftStimCtrl,2), lw=3, color=PHOTOSTIMCOLORS['laser_left'])
ax1.plot(0.3*np.array([-1,1])+4, 100*np.tile(meanRightStimCtrl,2), lw=3, color=PHOTOSTIMCOLORS['laser_right'])

(Z, pVal) = stats.ranksums(leftStimChange, controlLeftStimChange)
print 'p value between d1pi and control left hemi stim bias is {}'.format(pVal)
extraplots.significance_stars([1,2], 48, 2, starSize=10, gapFactor=0.12, color='0.5')
(Z, pVal) = stats.ranksums(rightStimChange, controlRightStimChange)
print 'p value between d1pi and control right hemi stim bias is {}'.format(pVal)
extraplots.significance_stars([3,4], 48, 2, starSize=10, gapFactor=0.12, color='0.5')

xlim = [0.5, 4.5]
ylim = [-50, 50]
plt.xlim(xlim)
plt.ylim(ylim)
xticks = [1,2,3,4]
xticklabels = ['Left\nstim', 'Control\nleft stim', 'Right\nstim', 'Control\nright stim']
plt.xticks(xticks, xticklabels, fontsize=fontSizeTicks)
plt.yticks(range(-40,60,10))
labelDis = 0.1
#plt.xlabel('Photostimulation', fontsize=fontSizeLabels) # labelpad=labelDis
plt.ylabel('Rightward bias (%)\nstim - control', fontsize=fontSizeLabels) # labelpad=labelDis
extraplots.boxoff(ax1)
ax1.spines['bottom'].set_visible(False)
[t.set_visible(False) for t in ax1.get_xticklines()]

#extraplots.significance_stars([1,2], 48, 2, starSize=10, gapFactor=0.12, color='0.5')



# -- Panel B: relationship between distance off from center and contralateral bias -- #
ax2 = plt.subplot(gs[0,1])
plt.hold(True)
extraplots.boxoff(ax2)
ax2.annotate('B', xy=(labelPosX[1],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


# -- Load data about behavior bias in photostim 2afc -- #
FIGNAME_behav = 'photostim_2afc'
dataDir_behav = os.path.join(settings.FIGURES_DATA_PATH, STUDY_NAME, FIGNAME_behav)

summaryFilename_behav = 'summary_photostim_percent_contra_choice_change.npz'
summaryFullPath_behav = os.path.join(dataDir_behav,summaryFilename_behav)
summary_behav = np.load(summaryFullPath_behav)

# These numbers are percent change in contralateral choice (stim - control) for each condition:
left014 = summary_behav['d1pi014leftHemiStim'][:maxSessionNum]
left015 = summary_behav['d1pi015leftHemiStim'][:maxSessionNum]
left016 = summary_behav['d1pi016leftHemiStim'][:maxSessionNum]
right014 = summary_behav['d1pi014rightHemiStim'][:maxSessionNum]
right015 = summary_behav['d1pi015rightHemiStim'][:maxSessionNum]
right016 = summary_behav['d1pi016rightHemiStim'][:maxSessionNum]

left014sessions = summary_behav['d1pi014leftHemiStimSessions']
left015sessions = summary_behav['d1pi015leftHemiStimSessions']
left016sessions = summary_behav['d1pi016leftHemiStimSessions']
right014sessions = summary_behav['d1pi014rightHemiStimSessions']
right015sessions = summary_behav['d1pi015rightHemiStimSessions']
right016sessions = summary_behav['d1pi016rightHemiStimSessions']


locationDict = {'0': {'d1pi014_left':left014,
                      'd1pi016_right':right016,
                      'd1pi015_left':left015,
                      'd1pi015_right':right015},
                '0.4': {'d1pi014_right':right014,
                        'd1pi016_left':left016}
                }

ax2.axhline(y=0, color='k', linestyle='-')
np.random.seed(7)
centerBias = []
borderBias = []
for key, valueDict in locationDict.items():
    for label, value in valueDict.items():
        if key == '0':
            centerBias.extend(value)
        elif key == '0.4':
            borderBias.extend(value)
        randOffset = 0.1*(np.random.rand(len(value))-0.5)
        if label.split('_')[1] == 'left':
            ax2.plot(float(key)+randOffset, 100*value, 'o', mec=PHOTOSTIMCOLORS['laser_left'], mfc='None')
        elif label.split('_')[1] == 'right':
            ax2.plot(float(key)+randOffset, 100*value, 'o', mec=PHOTOSTIMCOLORS['laser_right'], mfc='None')

meanCenterBias = np.mean(centerBias)
meanBorderBias = np.mean(borderBias)
ax2.plot(0.1*np.array([-0.5,0.5])+0, 100*np.tile(meanCenterBias,2), lw=2, color='grey')
ax2.plot(0.1*np.array([-0.5,0.5])+0.4, 100*np.tile(meanBorderBias,2), lw=2, color='grey')

ax2.set_xlim([-0.2, 0.6])
ax2.set_ylim([-15, 50])
ax2.set_xticks([0, 0.4])
ax2.set_xticklabels(['Center of\npost. striatum', 'Border of\npost. striatum'])
ax2.set_ylabel('Contralateral bias (%)\n stim - control', fontsize=fontSizeLabels)
#plt.show()
extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
extraplots.boxoff(ax2)
ax2.spines['bottom'].set_visible(False)
[t.set_visible(False) for t in ax2.get_xticklines()]

# -- Stats -- #
(Z, pVal) = stats.ranksums(centerBias, borderBias)
print 'p value between center and border stim sites is {}'.format(pVal)
extraplots.significance_stars([0,0.4], 48, 2, starSize=10, gapFactor=0.12, color='0.5')
if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)

