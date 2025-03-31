"""
Find the location of each recorded cell for a given subject (in CCF coordinates).

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

# Unfortunately, there is a bug with update in pandas v1.5.3 (though it works with v2.2.3)
# https://github.com/pandas-dev/pandas/issues/55509
celldb.update(celldbOneSubjectCoords) # The bug would change ints to floats
#celldb = celldb.combine_first(celldbOneSubjectCoords)

# BUG: for some reason 'recordingSiteName' is not being copied, so I added the line below:
#celldb.loc[celldbOneSubject.index].recordingSiteName = list(celldbOneSubject.recordingSiteName)
#celldb.recordingSiteName.loc[celldbOneSubject.index].update(celldbOneSubject.recordingSiteName)

celldatabase.save_hdf(celldb, dbCoordsFilename)

if 1:
    # -- Show all recording sites --
    aa = ha.AllenAverageCoronalAtlas()
    aa.add_points_from_db(celldb)
    aa.show_all_sites(nRows=2, areas=['AUDp','AUDv','AUDd','AUDpo'])

sys.exit()


# ---------------------------------------------------------------------------
#celldb[newColumns] = celldbOneSubjectCoords[newColumns]
#celldb['recordingSiteName'].fillna('', inplace=True)

#dbFilename = os.path.join(dbPath, f'celldb_coords.h5')
#celldatabase.save_hdf(celldb, dbFilename)
'''
subject = 'acid007'
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,'celldb_{}.h5'.format(studyparams.STUDY_NAME))
celldb = celldatabase.load_hdf(dbFilename)
celldbOneSubject = celldb[celldb.subject==subject]
'''



sys.exit()

'''
inforec = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')
basicdb = celldatabase.generate_cell_database(inforec, minimal=False)
celldb = ha.cell_locations(basicdb)
'''

#dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
#dbFilename = os.path.join(dbPath, f'celldb_coords_{subject}.h5')
#celldatabase.save_hdf(celldb, dbFilename)



subject = 'febe008' # 'test000' #'feat001'
inforec = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')

#dbPath = f'/tmp/celldb_{subject}.h5'

dbPath = os.path.join(settings.DATABASE_PATH, f'celldb_{subject}.h5')

sys.exit()

if 0:
    basicdb = celldatabase.generate_cell_database(inforec, minimal=False)
    celldatabase.save_hdf(basicdb, dbPath)
    sys.exit()

if 1:
    basicdb = celldatabase.load_hdf(dbPath)
    brainAreaDict = None #{'left_AudStr': 'LeftAstr', 'right_AudStr': 'RightAstr'}
    filterConditions = None
    celldb = ha.cell_locations(basicdb, filterConditions, brainAreaDict)
    celldatabase.save_hdf(celldb, dbPath)
    sys.exit()

celldb = celldatabase.load_hdf(dbPath)

#celldb = celldb[celldb.date=='2022-02-07']
# celldb.bestChannel.min(), celldb.bestChannel.max()


# -- Show all recording sites --
aa = ha.AllenAverageCoronalAtlas()
aa.add_points_from_db(celldb)
aa.show_all_sites(nRows=2, areas=['AUDp','AUDv','AUDd','AUDpo'])
#aa.show_all_sites()

figPath = os.path.join('/mnt/jarahubdata/reports/2022paspeech/feat_tracks/', f'{subject}_tracks_on_brain.png')
#plt.savefig('/tmp/feat006_tracks_on_brain.png', format='png')
plt.savefig(figPath, format='png')
