"""
Load and plot the waveform from all cells.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

dataDir = '/data/2018acsup_photoid/'
dfPath = os.path.join(dataDir,'dataframe_2018acsup_photoid_reduced.h5')
celldf = pd.read_hdf(dfPath)

# -- Create an array with the spike waveform for each cell --
spikeWaveforms = np.array(list(celldf.spikeShape))  # [nCells, nSamples]

spikeWaveformsNormalized = spikeWaveforms/abs(spikeWaveforms.min(axis=1)[:,np.newaxis])

plt.clf()
plt.plot(spikeWaveformsNormalized.T)
plt.show()
