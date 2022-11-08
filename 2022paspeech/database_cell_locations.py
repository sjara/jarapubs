"""
Find location of cells recorded by Jenny (feat001 mouse)
"""

import os
import sys
from jaratoolbox import celldatabase
from jaratoolbox import histologyanalysis as ha
import matplotlib.pyplot as plt
from jaratoolbox import settings
import studyparams
from importlib import reload
reload(ha)
reload(celldatabase)

print('input animal name (e.g. feat004)')
subject = input()
inforec = os.path.join(settings.INFOREC_PATH, f'{subject}_inforec.py')

#dbPath = f'/tmp/celldb_{subject}.h5'

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, f'celldb_{subject}.h5')

if 0:
    basicdb = celldatabase.generate_cell_database(inforec, minimal=False)
    celldatabase.save_hdf(basicdb, dbPath)
    sys.exit()

if 0:
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

figPath = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'feat_tracks' f'{subject}_tracks_on_brain.png')

#plt.savefig(figPath, format='png')
