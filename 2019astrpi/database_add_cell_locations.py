"""
Add cell locations to the database.
"""

import os
from jaratoolbox import celldatabase
from jaratoolbox import histologyanalysis as ha
from jaratoolbox import settings
import studyparams


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
origPath = os.path.join(figuresDataDir, 'astrpi_all_cells_no_stats.h5')
#origPath = '/data/figuresdata/2019astrpi/astrpi_all_cells_no_stats.h5'
dbPath = '/tmp/astrpi_all_cells_locations.h5'

basicDB = celldatabase.load_hdf(origPath)

brainAreaDict = {'left_AudStr': 'LeftAstr', 'right_AudStr': 'RightAstr'}
filterConditions = None
cellDB = ha.cell_locations(basicDB, filterConditions, brainAreaDict)

celldatabase.save_hdf(cellDB, dbPath)

# -- Show all recording sites --
if 0:
    aa = ha.AllenAverageCoronalAtlas()
    aa.add_points_from_db(cellDB)
    aa.show_all_sites()
