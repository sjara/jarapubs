"""
Estimate how many days passed between one session and the next.
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import datetime
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
import scipy.stats as stats 
import studyparams
import studyutils
import figparams

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

# -- Select only early or late days --
earlyDatesInds = [0,1]  # First two days
lateDatesInds = [2,3,4]   # Last two days
cellsEarly = np.zeros(len(celldb), dtype=bool)
cellsLate = np.zeros(len(celldb), dtype=bool)
datesEachSubject = []
for subject in celldb.subject.unique():
    subjectCells = celldb[celldb.subject==subject]
    datesThisSubject = subjectCells.date.unique()
    datesEachSubject.append(list(datesThisSubject))
    print(f'{subject}: {datesThisSubject}')
    for date in datesThisSubject:
        if date in datesThisSubject[earlyDatesInds]:
            cellsEarly = cellsEarly | ((celldb.subject==subject) & (celldb.date==date))
        elif date in datesThisSubject[lateDatesInds]:
            cellsLate = cellsLate | ((celldb.subject==subject) & (celldb.date==date))
        else:
            pass
celldbEarly = celldb.loc[cellsEarly].copy()
celldbLate = celldb.loc[cellsLate].copy()

diffAll = []
for inds, dates in enumerate(datesEachSubject):
    datesObj = [datetime.datetime.strptime(date, '%Y-%m-%d') for date in dates]
    differences = [(date2 - date1).days for date1, date2 in zip(datesObj[:-1], datesObj[1:])]
    diffAll.extend(differences)
    print(differences)

# -- Print median and quartiles --
print(f'Median: {np.median(diffAll)}')
print(f'Quartiles: {np.percentile(diffAll, [25, 75])}')
print(f'Min, Max: {np.min(diffAll)}, {np.max(diffAll)}')


