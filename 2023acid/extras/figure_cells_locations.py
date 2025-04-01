"""
Show where recorded neurons are located in the brain.
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
import scipy.stats as stats 
import studyparams
import studyutils
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(studyutils)
importlib.reload(studyparams)

def simplify_site_name(recordingSiteName):
    simplerName = recordingSiteName.split(',')[0]
    simplerName = simplerName.replace('Supplemental', 'Supp')
    if 'Ectorhinal' in simplerName:
        simplerName = 'Ectorhinal area'
    return simplerName

# -- Load data --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_coords.h5')
celldbCoords = celldatabase.load_hdf(dbFilename)
dbFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldb = celldatabase.load_hdf(dbFilename)
# NOTE: I had to create a new celldb file with the coordinates of each cell because
#       the original celldb file did not contain this information (bug in pandas).

recordingSiteName = celldbCoords.recordingSiteName.apply(simplify_site_name)
celldb.recordingSiteName = recordingSiteName

# Show counts
print(celldb.brainArea.groupby(recordingSiteName).count())
print('')

# -- Process data (copied from figure_overall_firing.py) --
celldbAll = celldb
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=5)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldbAll, steadyParams, maxChangeFactor)
includedCells = responsive & steady

celldbIncluded = celldbAll[includedCells]
# -- Calculate percentage of included cells per site --
countPerSite = celldbIncluded.brainArea.groupby(recordingSiteName).count()
print(countPerSite)
print('')
percentPerSite = countPerSite / countPerSite.sum() * 100
formatted = percentPerSite.map(lambda x: f"{x:.1f}%")
print(formatted)

auditoryAreas = ['auditory' in x for x in percentPerSite.index]
nonAuditoryAreas = [('auditory' not in x) and ('Temporal' not in x) for x in percentPerSite.index]
totalAud = percentPerSite[auditoryAreas].sum()
totalNonAud = percentPerSite[nonAuditoryAreas].sum()
totalTeA = percentPerSite['Temporal association areas']
total = totalAud + totalNonAud + totalTeA

print(f'\nPercentages:')
print(f'Auditory \t {totalAud:0.1f}%')
print(f'Non-auditory \t {totalNonAud:0.1f}%')
print(f'Temporal As. \t {totalTeA:0.1f}%')
print(f'Total \t {total}')

# 84 + 9.3 + 1.8 + 3.6 + 1.4 + 0.4
