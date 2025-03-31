"""
Generate and save database containing basic information, stats, and indices for each cell.
"""

import os
from jaratoolbox import settings
from jaratoolbox import celldatabase
import studyparams


if __name__ == "__main__":
    # -- Generate cell database (this function excludes bad cells by default, see docs) --
    celldb = celldatabase.generate_cell_database_from_subjects(studyparams.SUBJECTS)

    dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
    dbFilename = os.path.join(dbPath,'celldb_{}_basic.h5'.format(studyparams.STUDY_NAME))
    if os.path.isdir(dbPath):
        celldatabase.save_hdf(celldb, dbFilename)
        print('Saved database to {}'.format(dbFilename))
    else:
        print('{} does not exist. Please create this folder.'.format(dbPath))

