"""
Find the location of each recorded cell for a given subject (in CCF coordinates).

This script will update the celldb_STUDYNAME_coords.h5 database with cell locations
for one subject. If that database doesn't exist yet, it will create one based on
the basic database created by database_generation.py

You need to specify the subject to analyze, as in:
python database_cell_locations.py test000
"""

import os
import sys
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import histologyanalysis as ha  # This requires the Allen SDK
import matplotlib.pyplot as plt
from jaratoolbox import settings
import studyparams

from importlib import reload
reload(ha)
reload(celldatabase)

if len(sys.argv)==2:
    subject = sys.argv[1]
else:
    raise ValueError('You need to specify which subject to process.')

histoColumns = ['recordingSiteName', 'x_coord', 'y_coord', 'z_coord']
histoColumnsDefaultValue = ['', np.nan, np.nan, np.nan]

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbCoordsFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_coords.h5')
if os.path.isfile(dbCoordsFilename):
    celldb = celldatabase.load_hdf(dbCoordsFilename)
else:
    dbBasicFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_basic.h5')
    celldb = celldatabase.load_hdf(dbBasicFilename)
    histoColumnsDict = dict(zip(histoColumns, histoColumnsDefaultValue))
    celldb = celldb.assign(**histoColumnsDict)

celldbOneSubject = celldb[celldb.subject==subject]
celldbOneSubjectCoords = ha.cell_locations(celldbOneSubject)

# Unfortunately, there was a bug in update() in pandas v1.5.3 (though it works with v2.2.3)
# https://github.com/pandas-dev/pandas/issues/55509
# The bug would change ints to floats
celldb.update(celldbOneSubjectCoords)  # This works if the bug has been resolved
#celldb = celldb.combine_first(celldbOneSubjectCoords)  # Alternative if the bug exists

celldatabase.save_hdf(celldb, dbCoordsFilename)


# -- Show all recording sites --
aa = ha.AllenAverageCoronalAtlas()
aa.add_points_from_db(celldb)
aa.show_all_sites(nRows=1, areas=['AUDp','AUDv','AUDd','AUDpo'])

