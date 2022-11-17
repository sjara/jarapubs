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
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning.h5')
allSubjects = studyparams.EPHYS_MICE

# Check if fulldb (all subjects) exists, if doesn't exist, create it.
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
