"""
DESCRIBE HERE WHAT THE SCRIPT GENERATES.
"""

import os
import sys
import numpy as np
from jaratoolbox import settings
import studyparams

FIGNAME = 'figure_name'
figDataFile = 'file_containing_data_for_this_fig.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

# -- Optional --
if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)
    #print('Please create folder: {}'.format(figDataDir)); sys.exit()

figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

databaseName = 'THIS_STUDY_database.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
