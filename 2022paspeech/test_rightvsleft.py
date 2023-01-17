"""
Breaks results down by hemisphere
"""

import os
import sys
import numpy as np
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
from scipy import stats
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import studyparams
import figparams

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_leftVsRight' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [4, 6] # In inches
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')

audPColor = figparams.colors['audP']
audDColor = figparams.colors['audD']
audVColor = figparams.colors['audV']

selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])


ftCombos = np.array([selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax])
votCombos = np.array([selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax])
bestSelectivityIndexFt = ftCombos.max(0)
bestSelectivityIndexVot = votCombos.max(0)
isLeft = celldb.brainArea == 'AC_left'

FtSI_Left = np.full(len(celldb.brainArea), np.nan)
FtSI_Right = np.full(len(celldb.brainArea), np.nan)
FtSI_Left[isLeft] = bestSelectivityIndexFt[isLeft]
FtSI_Right[~isLeft] = bestSelectivityIndexFt[~isLeft]

VotSI_Left = np.full(len(celldb.brainArea), np.nan)
VotSI_Right = np.full(len(celldb.brainArea), np.nan)
VotSI_Left[isLeft] = bestSelectivityIndexVot[isLeft]
VotSI_Right[~isLeft] = bestSelectivityIndexVot[~isLeft]

bestFtSIbyArea_Left = []
bestFtSIbyArea_Right = []
bestVotSIbyArea_Left = []
bestVotSIbyArea_Right = []

hemispheres = np.unique(celldb.brainArea)

for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea_Left.append(FtSI_Left[recordingAreaName == thisArea])
    bestFtSIbyArea_Right.append(FtSI_Right[recordingAreaName == thisArea])
    bestVotSIbyArea_Left.append(VotSI_Left[recordingAreaName == thisArea])
    bestVotSIbyArea_Right.append(VotSI_Right[recordingAreaName == thisArea])

    bestFtSIbyArea_Left[indArea] = bestFtSIbyArea_Left[indArea][~np.isnan(bestFtSIbyArea_Left[indArea])]
    bestFtSIbyArea_Right[indArea] = bestFtSIbyArea_Right[indArea][~np.isnan(bestFtSIbyArea_Right[indArea])]
    bestVotSIbyArea_Left[indArea] = bestVotSIbyArea_Left[indArea][~np.isnan(bestVotSIbyArea_Left[indArea])]
    bestVotSIbyArea_Right[indArea] = bestVotSIbyArea_Right[indArea][~np.isnan(bestVotSIbyArea_Right[indArea])]

ustat, pValmannU_LvR_FtAudP = stats.mannwhitneyu(bestFtSIbyArea_Left[0], bestFtSIbyArea_Right[0])
ustat, pValmannU_LvR_FtAudD = stats.mannwhitneyu(bestFtSIbyArea_Left[1], bestFtSIbyArea_Right[1])
ustat, pValmannU_LvR_FtAudV = stats.mannwhitneyu(bestFtSIbyArea_Left[2], bestFtSIbyArea_Right[2])

ustat, pValmannU_LvR_VotAudP = stats.mannwhitneyu(bestVotSIbyArea_Left[0], bestVotSIbyArea_Right[0])
ustat, pValmannU_LvR_VotAudD = stats.mannwhitneyu(bestVotSIbyArea_Left[1], bestVotSIbyArea_Right[1])
ustat, pValmannU_LvR_VotAudV = stats.mannwhitneyu(bestVotSIbyArea_Left[2], bestVotSIbyArea_Right[2])
print('--Stats Summary--')
print('FT:')
print(f'AudP p = {pValmannU_LvR_FtAudP}')
print(f'AudD p = {pValmannU_LvR_FtAudD}')
print(f'AudV p = {pValmannU_LvR_FtAudV}')
print('VOT:')
print(f'AudP p = {pValmannU_LvR_VotAudP}')
print(f'AudD p = {pValmannU_LvR_VotAudD}')
print(f'AudV p = {pValmannU_LvR_VotAudV}')


plt.figure()
gsMain = gridspec.GridSpec(2,1)
gsMain.update(left=0.15, right=0.98, top=0.95, bottom=0.05, wspace=0.2, hspace=0.3)

ax1 = plt.subplot(gsMain[0,0])
labels = ['audD R', 'audD L', 'audP R', 'audP L', 'audV R', 'audV L']
plotData = [bestVotSIbyArea_Right[1], bestVotSIbyArea_Left[1], bestVotSIbyArea_Right[0], bestVotSIbyArea_Left[0], bestVotSIbyArea_Right[2], bestVotSIbyArea_Left[2]]
plt.boxplot(plotData, notch= True, labels = labels)
plt.ylabel('SI-VOT')
plt.title('VOT Selectivity x Area x hemisphere')

for indVars, theseVars in enumerate(plotData):
    y = plotData[indVars]
    x = np.random.normal(indVars+1, 0.04, size=len(y))
    if indVars < 2:
        colors = audDColor
    elif indVars > 3:
        colors = audVColor
    else:
        colors = audPColor
    plt.scatter(x,y, c = colors, edgecolors = 'r', alpha = 0.3)


ax2 = plt.subplot(gsMain[1,0])
plotData = [bestFtSIbyArea_Right[1], bestFtSIbyArea_Left[1], bestFtSIbyArea_Right[0], bestFtSIbyArea_Left[0], bestFtSIbyArea_Right[2], bestFtSIbyArea_Left[2]]
plt.boxplot(plotData, notch = True, labels = labels)
plt.ylabel('SI-FT')
plt.title('FT Selectivity x Area x hemisphere')
for indVars, theseVars in enumerate(plotData):
    y = plotData[indVars]
    x = np.random.normal(indVars+1, 0.04, size=len(y))
    if indVars < 2:
        colors = audDColor
    elif indVars > 3:
        colors = audVColor
    else:
        colors = audPColor
    plt.scatter(x,y, c = colors, edgecolors = 'r', alpha = 0.3)

plt.show()
extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
