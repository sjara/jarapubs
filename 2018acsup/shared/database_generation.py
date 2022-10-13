"""
Generate database (basic stats).

Based on ../database_generation_workflow.py 

This takes about 15 seconds.
"""

import sys
import pandas as pd
from jaratoolbox import celldatabase
sys.path.append('..')
import studyparams

from importlib import reload

dbPath = '/tmp/celldb_2018acsup_photoid_basic.h5'
#dbPath = '/data/figuresdata/2018acsup/shared/celldb_2018acsup_photoid_basic.h5'

celldbPV = celldatabase.generate_cell_database_from_subjects(studyparams.PV_CHR2_MICE)
celldbPV['mouseGenotype'] = 'PV::ChR2'
celldbSOM = celldatabase.generate_cell_database_from_subjects(studyparams.SOM_CHR2_MICE)
celldbSOM['mouseGenotype'] = 'SOM::ChR2'

celldbAll = pd.concat([celldbPV, celldbSOM], ignore_index=True)

# -- See database_photoidentification.py (line 221) and studyparams.py --
ISI_THRESHOLD = 0.02  # Maximum allowed fraction of ISI violations
SPIKE_QUALITY_THRESHOLD = 2.5  # Mininum quality

celldb = celldbAll.query('(isiViolations < @ISI_THRESHOLD) and '+
                         '(spikeShapeQuality > @SPIKE_QUALITY_THRESHOLD)')

# Fix cluster from float to int
celldb.cluster = celldb.cluster.astype(int)

'''
if 0:
    # -- The following line would create/overwrite .clu.modified files --
    from jaratoolbox import spikesorting
    celldb = spikesorting.rescue_clusters(celldb, isiThreshold=studyparams.ISI_THRESHOLD)
'''

celldatabase.save_hdf(celldb, dbPath)
