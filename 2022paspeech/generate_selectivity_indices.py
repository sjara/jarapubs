"""
Creates firing rate based selectivity indices for VOT and FT (maximum modulation by VOT for each FT) and compares these indices between cortical areas. Also checks following controls (comparing between areas): average of VOT SI between FTmin/FTmax, checking SI-VOT FTmin and SI-VOT FTmax separately, checking max(SI-VOT FTmin, SI-VOT FTmax) including all cells (not just cells that respond to speech sounds)
"""

import os
import sys
import numpy as np
from jaratoolbox import settings
import studyparams

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

# -- Optional --
if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)
    #print('Please create folder: {}'.format(figDataDir)); sys.exit()

figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')


selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])


ftCombos = ['selectivityIndexFT_VOTmin', 'selectivityIndexFT_VOTmax']
votCombos = ['selectivityIndexVOT_FTmin', 'selectivityIndexVOT_FTmax']
bestSelectivityIndexFt = celldb[[ftCombos[0], ftCombos[1]]].values.max(1)
bestSelectivitytIndexVot = celldb[[votCombos[0], votCombos[1]]].values.max(1)

bestFtSIbyArea = []
bestVotSIbyArea = []

for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea.append(bestSelectivityIndexFt[recordingAreaName == thisArea])
    bestVotSIbyArea.append(bestSelectivityIndexVot[recordingAreaName == thisArea])


np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityIndexFt = bestSelectivityIndexFt, bestSelectivitytIndexVot = bestSelectivitytIndexVot, audCtxAreas = audCtxAreas, bestFtSIbyArea = bestFtSIbyArea, bestVotSIbyArea = bestVotSIbyArea)
print('saved to ' f'{figDataFullPath}')
