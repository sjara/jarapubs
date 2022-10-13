"""
Test consistency of data. Check that npz files match cells in dataframe.

This will check mismatches produced when two sessions on the same day have the same pdepth.
"""

import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

dataDir = '/data/2018acsup_photoid/'
dfPath = os.path.join(dataDir,'dataframe_2018acsup_photoid_reduced.h5')
celldf = pd.read_hdf(dfPath)

plt.clf()
for indCell, (indRow, dfRow) in enumerate(celldf.iterrows()):

    dataFile = (f'{dfRow.subject}_{dfRow.date}_{dfRow.pdepth}um_' +
                f'g{dfRow.egroup}c{dfRow.cluster}.npz')

    # -- Load the spikes and laser timestamps data --
    tsData = np.load(os.path.join(dataDir,dataFile))

    # -- Make sure the waveform and timestamps are from the same cell --
    assert indRow==tsData['indRow']

    print(f'{indCell} [{indRow}] {dataFile}')
