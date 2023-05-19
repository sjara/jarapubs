"""
Test if there is a correlation between speech feature selectivity and selectivity to basic response properties
"""
import os
import sys
import studyparams
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

from importlib import reload
reload(extraplots)

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(databaseDir, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)

amOnsetSelective = np.array(celldb.amSelectivityPvalOnset < 0.05)
amSustainSelective = np.array(celldb.amSelectivityPvalSustain < 0.05)
amSelective = amOnsetSelective| amSustainSelective

toneSelective = np.array(celldb.toneSelectivityPval < 0.05)

audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)

## -- calculate selectivity indices
selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

## -- select best index b/w conditions
ftCombos = np.array([selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax])
votCombos = np.array([selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax])
bestSelectivityIndexFt = ftCombos.max(0)
bestSelectivityIndexVot = votCombos.max(0)

## -- run stats
bestFtSIbyArea = []
bestVotSIbyArea = []
speechResponsiveByArea = []
amSelectivebyArea = []
toneSelectivebyArea = []
for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea.append(bestSelectivityIndexFt[recordingAreaName == thisArea])
    bestVotSIbyArea.append(bestSelectivityIndexVot[recordingAreaName == thisArea])
    amSelectivebyArea.append(amSelective[recordingAreaName == thisArea])
    toneSelectivebyArea.append(toneSelective[recordingAreaName == thisArea])
    
kstat, pvalAmVot_AudP = stats.kruskal(bestVotSIbyArea[0][amSelectivebyArea[0] == 1],bestVotSIbyArea[0][amSelectivebyArea[0] == 0])
kstat, pvalAmVot_AudD = stats.kruskal(bestVotSIbyArea[1][amSelectivebyArea[1] == 1],bestVotSIbyArea[1][amSelectivebyArea[1] == 0])
kstat, pvalAmVot_AudV = stats.kruskal(bestVotSIbyArea[2][amSelectivebyArea[2] == 1],bestVotSIbyArea[2][amSelectivebyArea[2] == 0])

kstat, pvalToneFt_AudP = stats.kruskal(bestFtSIbyArea[0][toneSelectivebyArea[0]==1], bestFtSIbyArea[0][toneSelectivebyArea[0]==0])
kstat, pvalToneFt_AudD = stats.kruskal(bestFtSIbyArea[1][toneSelectivebyArea[1]==1], bestFtSIbyArea[1][toneSelectivebyArea[1]==0])
kstat, pvalToneFt_AudV = stats.kruskal(bestFtSIbyArea[2][toneSelectivebyArea[2]==1], bestFtSIbyArea[2][toneSelectivebyArea[2]==0])

