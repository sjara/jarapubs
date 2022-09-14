"""
Check if any feature shows spatial arrangement.
"""

import sys
sys.path.append('..')
import studyparams
import figparams
import os
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from jaratoolbox import extraplots
from jaratoolbox import spikesorting
from scipy import stats
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import studyutils

from importlib import reload
reload(figparams)
reload(studyutils)
reload(extraplots)

SAVE_FIGURE = 0
outputDir = '/tmp/'
#figFilename = 'fig5_am_responses' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [4.6, 3.8] # In inches
figSize = [7, 5.5] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
markerSize = 2  #figparams.rasterMarkerSize

#labelPosXtop = [0.01, 0.29, 0.52, 0.75]   # Horiz position for panel labels
#labelPosXbot = [0.01, 0.53, 0.75]   # Horiz position for panel labels
#labelPosY = [0.94, 0.45]    # Vert position for panel labels
labelPosXtop = [0.01, 0.5]   # Horiz position for panel labels
labelPosXbot = [0.01, 0.53, 0.75]   # Horiz position for panel labels
labelPosY = [0.97, 0.68, 0.32]    # Vert position for panel labels

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_am_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)
#extraDataPath = os.path.join(figuresDataDir, 'am_tuning.npz')
#amTuningData = np.load(extraDataPath)

pp = lambda x: f"['{x.subject}', '{x.date}', {x.pdepth}, {x.egroup}, {x.cluster}]"

cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']

def print_corr(feature, selected):
    corrSpearX, pValSpearX = stats.spearmanr(feature[selected], xc[selected])
    corrPearsX, pValPearsX = stats.pearsonr(feature[selected], xc[selected])
    corrSpearY, pValSpearY = stats.spearmanr(feature[selected], yc[selected])
    corrPearsY, pValPearsY = stats.pearsonr(feature[selected], yc[selected])
    corrSpearZ, pValSpearZ = stats.spearmanr(feature[selected], zc[selected])
    corrPearsZ, pValPearsZ = stats.pearsonr(feature[selected], zc[selected])

    print(f'N = {len(feature[selected])}')
    print(f'Spearman X r={corrSpearX:0.2f} (p={pValSpearX:0.4f})')
    print(f'Pearson  X r={corrPearsX:0.2f} (p={pValPearsX:0.4f})')
    print('')
    print(f'Spearman Y r={corrSpearY:0.2f} (p={pValSpearY:0.4f})')
    print(f'Pearson  Y r={corrPearsY:0.2f} (p={pValPearsY:0.4f})')
    print('')
    print(f'Spearman Z r={corrSpearZ:0.2f} (p={pValSpearZ:0.4f})')
    print(f'Pearson  Z r={corrPearsZ:0.2f} (p={pValPearsZ:0.4f})')
    print('')


# -- Tone response index --
restrictND1 = 0  # True: use only ND1 from tetrodes with D1
toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1 = \
    studyutils.select_cells(celldb, restrictND1=restrictND1)

xc, yc, zc = celldb.x_coord.to_numpy(), celldb.y_coord, celldb.z_coord
midP = 227 # 227
xc[xc>midP] = midP - (xc[xc>midP]-midP)  # Flip right hemisphere
#xc[xc>midP] = np.nan  # Exclude the other hemisphere
#xc[(zc<255)|(zc>274)] = np.nan  # Exclude anterior or posterior
#xc[zc>274] = np.nan  # Exclude anterior or posterior
#xc[celldb.subject!='d1pi046'] = np.nan  # Exclude anterior or posterior

#subset = (indD1 | indND1); print('================ Using both D1 & non-D1 =================')
subset = indD1; print('================ Using D1 =================')
#subset = indND1; print('================ Using non-D1 =================')

validCF = ~np.isnan(celldb.toneCharactFreq)
goodFit = celldb.toneGaussianRsquare > 0.01
intThresh = celldb.toneIntensityThreshold 

# -- Threshold --
selAll = toneResponsive & validCF & ~np.isnan(xc)
selAll = toneResponsive & validCF & ~np.isnan(xc) & subset
print('-- Threshold --')
print_corr(celldb.toneIntensityThreshold, selAll)


# -- Best AM rate --
selAll = amSustainResponsive & ~np.isnan(xc)
selAll = amSustainResponsive & ~np.isnan(xc) & subset
print('-- Best AM rate --')
print_corr(celldb.amIndexBestSustain, selAll)


# -- Characteristic freq --
selAll = toneResponsive & validCF  & ~np.isnan(xc)
selAll = toneResponsive & validCF  & ~np.isnan(xc) & subset
#selAll = toneResponsive & goodFit & validCF & ~np.isnan(xc) # & (intThresh<60)
#selAll = toneResponsive & goodFit & validCF & ~np.isnan(xc) & indD1
#selD1 = indD1 & toneResponsive & goodFit & validCF
#selND1 = indND1 & toneResponsive & goodFit & validCF

logCharactFreq = np.log2(celldb.toneCharactFreq)
#logCharactFreq = (celldb.toneIntensityThreshold)
#logCharactFreq = (celldb.toneLatencySoundDetector)
#logCharactFreq = celldb.toneIndexBest
#logCharactFreq = celldb.amIndexBestSustain

print('-- Characteristic freq --')
print_corr(logCharactFreq, selAll)


# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(1, 3)
gsMain.update(left=0.09, right=0.96, top=0.9, bottom=0.1, hspace=0.35)


np.random.seed(1)
jitter = 0.2*np.random.rand(np.sum(selAll))
ylims = None #[11,14]
plt.subplot(gsMain[0])
plt.plot(xc[selAll], logCharactFreq[selAll]+jitter, '.')
#plt.gca().set_yscale('log')
#plt.ylim(ylims)
plt.xlim([65,125]) 
plt.title('X coord')
plt.subplot(gsMain[1])
plt.plot(yc[selAll], logCharactFreq[selAll]+jitter, '.')
plt.ylim(ylims)
plt.title('Y coord')
plt.subplot(gsMain[2])
plt.plot(zc[selAll], logCharactFreq[selAll]+jitter, '.')
plt.ylim(ylims)
plt.title('Z coord')
plt.show()

sys.exit()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
