'''
Figure of firing rate based selectivity indices for VOT and FT
'''

import os
import sys
import studyparams
import figparams
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
from scipy import signal
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from jaratoolbox import extraplots
import studyutils
from importlib import reload
reload(studyutils)
reload(studyparams)

SAVE_FIGURE = 0
outputDir = 'C:\\Users\\jenny\\tmp'
FIGNAME1 = 'figure_speech_selectivity_vot_allcells'
FIGNAME2 = 'figure_speech_selectivity_ft_allcells'
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches
#figDataFile = 'file_containing_data_for_this_fig.npz'
#figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
#figDataFullPath = os.path.join(figDataDir, figDataFile)

# -- Load data --
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir, figDataFile)
figData = np.load(figDataFullPath)
#databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning_allcells.h5')
#celldb = celldatabase.load_hdf(dbPath)
#audCtxAreas = ['Primary auditory area','Posterior auditory area', 'Ventral auditory area']


PANELS = [0] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.36, 0.7]   # Horiz position for panel labels
labelPosY = [0.9, 0.48]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
audPColor = figparams.colors['audP']
audDColor = figparams.colors['audD']
audVColor = figparams.colors['audV']

# -- Processed data --
kstat, pValKruskalBestVOT = stats.kruskal(bestVotSIbyArea[0], bestVotSIbyArea[1], bestVotSIbyArea[2], nan_policy = 'omit')
kstat, pValKruskalBestFT = stats.kruskal(bestFtSIbyArea[0], bestFtSIbyArea[1], bestFtSIbyArea[2], nan_policy = 'omit')

# if group difference, test individual compariosns:
if pValKruskalBestVOT < 0.05:
    ustat, pValmannU_votAudPvsAudD = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[1])
    ustat, pValmannU_votAudPvsAudV = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[2])
    ustat, pValmannU_votAudDvsAudV = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[2])

if pValKruskalBestFT < 0.05:
    ustat, pValmannU_ftAudPvsAudD = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[1])
    ustat, pValmannU_ftAudPvsAudV = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudDvsAudV = stats.mannwhitneyu(bestGtSIbyArea[1], bestFtSIbyArea[2])


# -- Plot VOT results --
fig1 = plt.gcf()
fig1.clf()
fig1.set_facecolor('w')

gsMain = gridspec.GridSpec(3,1)
gsMain.update(left=0.15, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)
#bins = np.arange(0,1,0.025)
#yMax = 14
bins = np.arange(0,1,0.05)
yMax = 18


ax1 = plt.subplot(gsMain[0, 0])
plt.hist(bestVotSIbyArea[1], bins = bins, color = audDColor)
plt.ylim([0, yMax])
audD_votmedian = np.nanmedian(bestVotSIbyArea[1])
plt.plot([audD_votmedian, audD_votmedian], [0,20], color = 'k', ls = '--')

ax2 = plt.subplot(gsMain[1, 0], sharex = ax1)
plt.hist(bestVotSIbyArea[0], bins = bins, color = audPColor)
plt.ylim([0, yMax])
audP_votmedian = np.nanmedian(bestVotSIbyArea[0])
plt.plot([audP_votmedian, audP_votmedian], [0,20], color = 'k', ls = '--')

ax3 = plt.subplot(gsMain[2, 0], sharex = ax1)
plt.hist(bestVotSIbyArea[2], bins = bins, color = audVColor)
plt.ylim([0, yMax])
audV_votmedian = np.nanmedian(bestVotSIbyArea[2])
plt.plot([audV_votmedian, audV_votmedian], [0,20], color = 'k', ls = '--')
#ax1.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

plt.show()
if SAVE_FIGURE:
    extraplots.save_figure(FIGNAME1, figFormat, figSize, outputDir)

# -- Plot FT results --
plt.figure()
fig2 = plt.gcf()
fig2.clf()
fig2.set_facecolor('w')

gsMain2 = gridspec.GridSpec(3,1)
gsMain2.update(left=0.15, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)
#bins = np.arange(0,1,0.025)
#yMax = 14
bins = np.arange(0,1,0.05)
yMax = 18


ax1 = plt.subplot(gsMain2[0, 0])
plt.hist(bestFtSIbyArea[1], bins = bins, color = audDColor)
plt.ylim([0, yMax])
audD_ftmedian = np.nanmedian(bestFtSIbyArea[1])
plt.plot([audD_ftmedian, audD_ftmedian], [0,20], color = 'k', ls = '--')

ax2 = plt.subplot(gsMain2[1, 0], sharex = ax1)
plt.hist(bestFtSIbyArea[0], bins = bins, color = audPColor)
plt.ylim([0, yMax])
audP_ftmedian = np.nanmedian(bestFtSIbyArea[0])
plt.plot([audP_ftmedian, audP_ftmedian], [0,20], color = 'k', ls = '--')

ax3 = plt.subplot(gsMain2[2, 0], sharex = ax1)
plt.hist(bestFtSIbyArea[2], bins = bins, color = audVColor)
plt.ylim([0, yMax])
audV_ftmedian = np.nanmedian(bestFtSIbyArea[2])
plt.plot([audV_ftmedian, audV_ftmedian], [0,20], color = 'k', ls = '--')

plt.show()

if SAVE_FIGURE:
    plt.figure(fig2)
    extraplots.save_figure(FIGNAME2, figFormat, figSize, outputDir)
