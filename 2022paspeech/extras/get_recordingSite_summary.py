import os
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
import csv

import sys
sys.path.append('..')
import studyparams

subject = 'feat010'
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, f'celldb_{subject}.h5')

celldb = celldatabase.load_hdf(dbPath)
sessionName = f'{subject} {celldb.date[1]} {celldb.maxDepth[1]}um'
summaryPath = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'cell_reports', 'feat_recording_sites.csv')

csvFile = open(summaryPath, 'a')
fieldNames = ['subject','date','maxDepth','recordingSiteName','x','y','z']
writer = csv.DictWriter(csvFile, fieldnames = fieldNames)
#writer.writeheader()
for indRow, dbRow in celldb.iterrows():
    writer.writerow({'subject':f'{dbRow.subject}', 'date':f'{dbRow.date}', 'maxDepth':f'{dbRow.maxDepth}', 'recordingSiteName':f'{dbRow.recordingSiteName}', 'x':f'{dbRow.x_coord}', 'y':f'{dbRow.y_coord}', 'z':f'{dbRow.z_coord}'})

csvFile.close()
