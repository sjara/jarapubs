"""
Select only cells with sound sessions (tones or AM).
"""

import sys
import os
import studyparams
from jaratoolbox import celldatabase
from jaratoolbox import settings

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
origdbPath = os.path.join(figuresDataDir, 'astrpi_all_cells_locations.h5')

newdbPath = '/tmp/astrpi_cells_with_soundsessions.h5'


origdb = celldatabase.load_hdf(origdbPath)
newdb = origdb.copy()

# Selects cells with a laserpulse session for D1 vs. nD1 determination
newdb = newdb[newdb['sessionType'].apply(lambda s: 'laserpulse' in s)]    

# Selects cells with either tuning curve or AM session for sound comparison
newdb = newdb[newdb['sessionType'].apply(lambda s: ('tuningCurve' in s) or ('am' in s))]

celldatabase.save_hdf(newdb, newdbPath)
