"""
Create a smaller database by excluding bad cells according to "manually_removed" file.
"""

import os
import sys
import studyparams
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings

toExclude = np.loadtxt(os.path.join('./', 'sj_manually_removed.txt'), dtype=int)

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
origdbPath = os.path.join(figuresDataDir, 'astrpi_cells_with_soundsessions.h5')
origdb = celldatabase.load_hdf(origdbPath)

newdbPath = '/tmp/astrpi_included_cells.h5'

newdb = origdb.drop(toExclude)
celldatabase.save_hdf(newdb, newdbPath)

