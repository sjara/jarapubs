import os
import sys
import studyparams
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
from scipy import signal

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

import studyutils

from importlib import reload
reload(studyutils)

# Load the data

databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning.h5')
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning_combineAudDAudPo.h5')
allSubjects = studyparams.EPHYS_MICE

# Try loading fulldb (all subjects), if doesn't exist, create it.
try:
    celldb = celldatabase.load_hdf(dbPath)
except:
    for indMouse, thisMouse in enumerate(allSubjects):
        subject = thisMouse
        dbPath = os.path.join(databaseDir, f'{subject}_paspeech_speech_pval.h5')
        if indMouse == 0:
            celldb = celldatabase.load_hdf(dbPath)
        else:
            onedb = celldatabase.load_hdf(dbPath)
            celldb =  celldb.append(onedb, ignore_index=True)
    recordingAreaName = celldb.recordingSiteName
    recordingAreaName = recordingAreaName.str.replace(r'[,/].+', '', regex = True)

    celldb['recordingAreaName'] = recordingAreaName
    newdbPath = os.path.join(databaseDir, 'fulldb_speech_tuning.h5')
    celldatabase.save_hdf(celldb, newdbPath)




'''
#audCtxAreas = ['Primary auditory area','Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Posterior auditory area', 'Ventral auditory area']
featureAlpha = 0.05/4
VOTselective_FTmax = celldb.votSelectivityFtMaxPvalOnset < featureAlpha
VOTselective_FTmin = celldb.votSelectivityFtMinPvalOnset < featureAlpha
FTselective_VOTmax = celldb.ftSelectivityVotMaxPvalOnset < featureAlpha
FTselective_VOTmin = celldb.ftSelectivityVotMinPvalOnset < featureAlpha
VOTnoFT = (VOTselective_FTmax | VOTselective_FTmin) & (~FTselective_VOTmax & ~FTselective_VOTmin)
FTnoVOT = (FTselective_VOTmax | FTselective_VOTmin) & (~VOTselective_FTmax & ~VOTselective_FTmin)
mixedSelective = (FTselective_VOTmax | FTselective_VOTmin) & (VOTselective_FTmax | VOTselective_FTmin)
notSelective = ~FTselective_VOTmax & ~FTselective_VOTmin & ~VOTselective_FTmax & ~VOTselective_FTmin
#bothFTandVOTselective = np.sum((VOTselective_FTmax == 1 | VOTselective_FTmin == 1) & (FTselective_VOTmax==1 | FTselective_VOTmin == 1)) #selective to VOT/FT for either irrelevant feature dimension. Should this alpha be 0.05 or does this need to be corrected for multiple comparisons (I think corrected)?

#selectivityByArea = zeros([len(audCtxAreas), 5])
mixedSelectiveByArea = np.empty([len(allSubjects), len(audCtxAreas), 5])

for indMouse, thisMouse in enumerate(allSubjects):
    for indArea, thisArea in enumerate(audCtxAreas):
        cellsThisMouse = celldb.subject == thisMouse
        cellsThisArea = celldb.recordingAreaName == thisArea
        numCells = np.sum(cellsThisMouse & cellsThisArea)
        #numVOTselectiveFtMax = np.sum(VOTselective_FTmax & (celldb.recordingAreaName == thisArea))
        #numVOTselectiveFtMin = np.sum(VOTselective_FTmin & (celldb.recordingAreaName == thisArea))
        #numFTselectiveVotMax = np.sum(FTselective_VOTmax & (celldb.recordingAreaName == thisArea))
        #numFTselectiveVotMin = np.sum(FTselective_VOTmin & (celldb.recordingAreaName == thisArea))
        #numBoth = np.sum(bothFeats[celldb.recordingAreaName == thisArea])
        #numFTOnly = np.sum(FTnoVOT[celldb.recordingAreaName == thisArea])
        #numVOTOnly = np.sum(VOTnoFT[celldb.recordingAreaName == thisArea])
        #numNeither = np.sum((FTnoVOT[celldb.recordingAreaName == thisArea]==0) & (VOTnoFT[celldb.recordingAreaName == thisArea]==0) & (bothFeats[celldb.recordingAreaName == thisArea] == 0) )
        #numAny =  np.sum(VOTselective_FTmin[celldb.recordingAreaName == thisArea] | VOTselective_FTmax[celldb.recordingAreaName == thisArea] | FTselective_VOTmin[celldb.recordingAreaName == thisArea] | FTselective_VOTmax[celldb.recordingAreaName == thisArea])
        numVOTOnly = np.sum(VOTnoFT[cellsThisArea & cellsThisMouse])
        numFTOnly = np.sum(FTnoVOT[cellsThisArea & cellsThisMouse])
        numBoth = np.sum(mixedSelective[cellsThisArea & cellsThisMouse])
        numNeither = np.sum(notSelective[cellsThisArea & cellsThisMouse])

        mixedSelectiveByArea[indMouse, indArea, 0] = numBoth
        mixedSelectiveByArea[indMouse, indArea, 1] = numFTOnly
        mixedSelectiveByArea[indMouse, indArea, 2] = numVOTOnly
        mixedSelectiveByArea[indMouse, indArea, 3] = numNeither
        mixedSelectiveByArea[indMouse, indArea, 4] = numCells
        #mixedSelectiveByArea[indArea, 5] = numAny
'''
    '''
    selectivityByArea[indArea, 0] = numVOTselectiveFtMax
    selectivityByArea[indArea, 1] = numVOTselectiveFtMin
    selectivityByArea[indArea, 2] = numFTselectiveVotMax
    selectivityByArea[indArea, 3] = numFTselectiveVotMin
    selectivityByArea[indArea, 4] = numCells
    '''
'''
expectedSelectiveByArea = zeros([len(audCtxAreas), 4])
sumCategories = np.sum(mixedSelectiveByArea, axis = 0)
for indArea, thisArea in enumerate(audCtxAreas):
    proportionTotalCells = mixedSelectiveByArea[indArea,4]/sumCategories[4]
    expectBoth =  sumCategories[0]*proportionTotalCells
    expectFTonly = sumCategories[1]*proportionTotalCells
    expectVOTonly = sumCategories[2]*proportionTotalCells
    expectNeither = sumCategories[3]*proportionTotalCells

    expectedSelectiveByArea[indArea, 0] = expectBoth
    expectedSelectiveByArea[indArea, 1] = expectFTonly
    expectedSelectiveByArea[indArea, 2] = expectVOTonly
    expectedSelectiveByArea[indArea, 3] = expectNeither






sessionSummary = zeros([2,len(eachSession)])
for indSession, thisSession in enumerate(eachSession):
    totalNumCells = np.sum((celldb.date == thisSession))
    numAudPCells = np.sum((celldb.date == thisSession) & (celldb.recordingAreaName == audCtxAreas[0]))
    numAudPoCells = np.sum((celldb.date == thisSession) & (celldb.recordingAreaName == audCtxAreas[1]))
    numAudDCells = np.sum((celldb.date == thisSession) & (celldb.recordingAreaName == audCtxAreas[2]))
    numAudVCells = np.sum((celldb.date == thisSession) & (celldb.recordingAreaName == audCtxAreas[2]))
    numAudCells = numAudPCells + numAudPoCells + numAudDCells + numAudVCells
    sessionSummary[0,indSession] = totalNumCells
    sessionSummary[1,indSession] = numAudCells
'''
