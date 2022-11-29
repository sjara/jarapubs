'''
Figure of firing rate based selectivity indices for VOT and FT
'''
import os
import sys
import studyparams
import figparams
import numpy as np
import pandas as pd
from jaratoolbox import settings
from jaratoolbox import extraplots
from scipy import stats
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from importlib import reload
reload(studyparams)

SAVE_FIGURE = 1
outputDir = 'C:\\Users\\jenny\\tmp'
FIGNAME = 'selectivityIndices'
FIGNAME1 = 'figure_speech_selectivity_vot_allcells'
FIGNAME2 = 'figure_speech_selectivity_ft_allcells'
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [9.5, 8.5] # In inches

# -- Load data --
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)
#databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning_allcells.h5')
#audCtxAreas = ['Primary auditory area','Posterior auditory area', 'Ventral auditory area']


PANELS = [0] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.48, 0.7]   # Horiz position for panel labels
labelPosY = [0.9, 0.48]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
audPColor = figparams.colors['audP']
audDColor = figparams.colors['audD']
audVColor = figparams.colors['audV']


## NEED TO UPDATE THIS TO WORK WITH NPZ FILE!
bestFtSIbyArea = []
bestVotSIbyArea = []
speechResponsiveByArea = []

for indArea, thisArea in enumerate(figData['audCtxAreas']):
    bestFtSIbyArea.append(figData['bestSelectivityIndexFt'][figData['recordingAreaName'] == thisArea])
    bestVotSIbyArea.append(figData['bestSelectivityIndexVot'][figData['recordingAreaName'] == thisArea])
    speechResponsiveByArea.append(figData['speechResponsive'][figData['recordingAreaName'] == thisArea])

# -- if group difference, test individual comparisons:
## all cells:
if figData['pValKruskalBestVOT'] < 0.05:
    ustat, pValmannU_votAudPvsAudD = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[1])
    ustat, pValmannU_votAudPvsAudV = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[2])
    ustat, pValmannU_votAudDvsAudV = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[2])

if figData['pValKruskalBestFT'] < 0.05:
    ustat, pValmannU_ftAudPvsAudD = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[1])
    ustat, pValmannU_ftAudPvsAudV = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudDvsAudV = stats.mannwhitneyu(bestGtSIbyArea[1], bestFtSIbyArea[2])

## responsive cells
if figData['pValKruskalBestVOT'] < 0.05:
    ustat, pValmannU_votAudPvsAudD = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_votAudPvsAudV = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_votAudDvsAudV = stats.mannwhitneyu(bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2]])

if figData['pValKruskalBestFT'] < 0.05:
    ustat, pValmannU_ftAudPvsAudD = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_ftAudPvsAudV = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_ftAudDvsAudV = stats.mannwhitneyu(bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2]])


#fig1 = plt.gcf()
#fig1.clf()
fig1 = plt.figure()
gsMain = gridspec.GridSpec(3,2)
gsMain.update(left=0.08, right=0.98, top=0.95, bottom=0.1, wspace=0.4, hspace=0.4)

bins = np.arange(0,1,0.05)
yMax = 25

# -- Plot VOT results --
ax1 = plt.subplot(gsMain[0, 0])
plt.hist(bestVotSIbyArea[1][speechResponsiveByArea[1]], bins = bins, color = audDColor)
plt.ylim([0, yMax])
audD_votmedian = np.nanmedian(bestVotSIbyArea[1][speechResponsiveByArea[1]])
plt.plot([audD_votmedian, audD_votmedian], [0,yMax], color = 'k', ls = '--')
plt.title('VOT Selectivity', fontsize=fontSizeLabels, fontweight='bold')
ax1.annotate(f'med = {np.round(audD_votmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

ax2 = plt.subplot(gsMain[1, 0], sharex = ax1)
plt.hist(bestVotSIbyArea[0][speechResponsiveByArea[0]], bins = bins, color = audPColor)
plt.ylim([0, yMax])
audP_votmedian = np.nanmedian(bestVotSIbyArea[0][speechResponsiveByArea[0]])
plt.plot([audP_votmedian, audP_votmedian], [0,yMax], color = 'k', ls = '--')
plt.ylabel('Cell Count', fontsize = fontSizeTicks)
ax2.annotate(f'med = {np.round(audP_votmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

ax3 = plt.subplot(gsMain[2, 0], sharex = ax1)
plt.hist(bestVotSIbyArea[2][speechResponsiveByArea[2]], bins = bins, color = audVColor)
plt.ylim([0, yMax])
audV_votmedian = np.nanmedian(bestVotSIbyArea[2][speechResponsiveByArea[2]])
plt.plot([audV_votmedian, audV_votmedian], [0,yMax], color = 'k', ls = '--')
plt.xlabel('VOT Selectivity Index', fontsize=fontSizeTicks)
ax3.annotate(f'med = {np.round(audV_votmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

# -- Plot FT results --
ax4 = plt.subplot(gsMain[0, 1])
plt.hist(bestFtSIbyArea[1][speechResponsiveByArea[1]], bins = bins, color = audDColor)
plt.ylim([0, yMax])
audD_ftmedian = np.nanmedian(bestFtSIbyArea[1][speechResponsiveByArea[1]])
plt.plot([audD_ftmedian, audD_ftmedian], [0,yMax], color = 'k', ls = '--')
plt.title('FT Selectivity', fontsize=fontSizeLabels, fontweight='bold')
ax4.annotate(f'med = {np.round(audD_ftmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

ax5 = plt.subplot(gsMain[1, 1], sharex = ax1)
plt.hist(bestFtSIbyArea[0][speechResponsiveByArea[0]], bins = bins, color = audPColor)
plt.ylim([0, yMax])
audP_ftmedian = np.nanmedian(bestFtSIbyArea[0][speechResponsiveByArea[0]])
plt.plot([audP_ftmedian, audP_ftmedian], [0,yMax], color = 'k', ls = '--')
ax5.annotate(f'med = {np.round(audP_ftmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

ax6 = plt.subplot(gsMain[2, 1], sharex = ax1)
plt.hist(bestFtSIbyArea[2][speechResponsiveByArea[2]], bins = bins, color = audVColor)
plt.ylim([0, yMax])
audV_ftmedian = np.nanmedian(bestFtSIbyArea[2][speechResponsiveByArea[2]])
plt.plot([audV_ftmedian, audV_ftmedian], [0,yMax], color = 'k', ls = '--')
plt.xlabel('FT Selectivity Index', fontsize=fontSizeTicks)
ax6.annotate(f'med = {np.round(audV_ftmedian,2)}', xy=(0.6, 15), xycoords = 'data', fontsize = fontSizeTicks)

ax1.annotate('AudD', xy=(0.48, 0.82), xycoords='figure fraction', fontsize=fontSizeLabels, fontweight='bold')
ax2.annotate('AudP', xy=(0.48, 0.51), xycoords='figure fraction', fontsize=fontSizeLabels, fontweight='bold')
ax3.annotate('AudV', xy=(0.48, 0.20), xycoords='figure fraction', fontsize=fontSizeLabels, fontweight='bold')

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(FIGNAME, figFormat, figSize, outputDir)
