"""
Plot histogram of average spike amplitudes and trace standard deviations.
"""

import os
import numpy as np
from matplotlib import pyplot as plt
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import extraplots
import studyparams
import figparams
import importlib
importlib.reload(figparams)
importlib.reload(celldatabase)

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_spikeSNR.h5')
celldb = celldatabase.load_hdf(dbFilename)

colorAmp = figparams.colorVoltageTrace

xLim = 420

# -- Plot results --
plt.clf()
#plt.hist(celldb['meanAmplitude'], bins=40)
plt.hist(celldb.stDevTrace[celldb.rawDataExists], 11, fc='0.7', ec='w', alpha=1);
plt.hist(celldb['meanAmplitude'], bins=40, fc=colorAmp, ec='w', alpha=1)
#plt.xlabel('Mean amplitude or StDev of trace (uV)')
plt.xlabel('Voltage (μV)', fontsize=figparams.fontSizeLabels)
plt.ylabel('N neurons', fontsize=figparams.fontSizeLabels)
plt.ylim(0, 80)
plt.xlim(0, xLim)
extraplots.boxoff(plt.gca())
#plt.legend(['Voltage trace Std Dev', 'Mean spike amplitude'])
plt.legend(['Voltage trace σ', 'Spike amplitude'], fontsize=figparams.fontSizeLabels)
plt.show()

print(f'N neurons above xlim ({xLim}μV):', np.sum(celldb['meanAmplitude']>xLim))

extraplots.save_figure('plots_spike_amplitude_histogram', 'svg', [6, 2], '/tmp/')
