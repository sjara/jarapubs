"""
Show response of pStr neuron during AC terminal inactivation.
"""

import os
import sys
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from jaratoolbox import spikesanalysis
from jaratoolbox import extraplots
from jaratoolbox import spikesorting
from jaratoolbox import colorpalette as cp

STUDY_NAME = '2021terminactiv'
dbFilePath = os.path.join(settings.DATABASE_PATH, STUDY_NAME, 'celldb.h5')
celldb = celldatabase.load_hdf(dbFilePath)


outputDir = f'/data/reports/{STUDY_NAME}/cellreports/'

#laserColors = 20*['g','lightgreen']
laserColor = cp.TangoPalette['Chameleon3']  

cellDict = {'subject' : 'arch004',
            'date' : '2021-05-26',
            'depth' : 3000,
            'tetrode' : 6,
            'cluster' : 4}

indRow, dbRow = celldatabase.find_cell(celldb, **cellDict)
oneCell = ephyscore.Cell(dbRow)

# -- Prepare figure --
plt.clf()

gs = gridspec.GridSpec(2, 2, width_ratios=[1.5,1])
gs.update(left=0.07, right=0.98, top=0.97, bottom=0.17, wspace=0.4, hspace=0.2)

axFreqTuningBoth = plt.subplot(gs[0:2, 0])
axFreqTuningNoLaser = plt.subplot(gs[0, 1])
axFreqTuningLaser = plt.subplot(gs[1, 1])

markerSize = 5
fontSizeLabel = 14
fontSizeTicks = 13

# ------------------------ Frequency tuning ----------------------------
ephysData, bdata = oneCell.load('laserTuningCurve')
spikeTimes = ephysData['spikeTimes']
eventOnsetTimes = ephysData['events']['stimOn']
nTrials = len(bdata['currentFreq'])
eventOnsetTimes = eventOnsetTimes[0:nTrials]
timeRange = [-0.1, 0.3]  # In seconds

# -- Analyze responses --
(spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
    spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

paramEachTrial = bdata['currentFreq']
possibleParams = np.unique(bdata['currentFreq'])
laserTrial = bdata['laserTrial']
possibleLaser = [0, 1]
trialsEachCond = behavioranalysis.find_trials_each_combination(paramEachTrial,
                                                               possibleParams,
                                                               laserTrial, possibleLaser)
# -- Calculate tuning curve --
responsePeriod = [0, 0.1]
responsePeriodDuration = responsePeriod[-1]-responsePeriod[0]
laserTrialBool = laserTrial.astype(bool)
spkT2C = spikesanalysis.spiketimes_to_spikecounts
spikeCountMat = spkT2C(spikeTimesFromEventOnset,
                       indexLimitsEachTrial,
                       responsePeriod)
nParams = len(possibleParams)
nLaser = len(possibleLaser)
meanSpikeCountEachCond = np.empty((nParams, nLaser))
stdSpikeCountEachCond = np.empty((nParams, nLaser))
for indLaser,thisLaser in enumerate(possibleLaser):
    for indParam,thisParam in enumerate(possibleParams):
        thisSpkCount = spikeCountMat[trialsEachCond[:,indParam,indLaser]]
        meanSpikeCountEachCond[indParam,indLaser] = np.mean(thisSpkCount)
        stdSpikeCountEachCond[indParam,indLaser] = np.std(thisSpkCount)

meanFiringRateEachCond = meanSpikeCountEachCond/responsePeriodDuration
stdFiringRateEachCond = stdSpikeCountEachCond/responsePeriodDuration

#freqLabels = [int(possibleParams[0])] + (len(possibleParams)-2)*[None] + [int(possibleParams[-1])] 
freqLabels = [int(possibleParams[0]/1000)] + (len(possibleParams)-2)*[None] + [int(possibleParams[-1]/1000)] 

# -------------------------------------------------------------------------------
plt.sca(axFreqTuningNoLaser)
#pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange,
#                                               trialsEachCond[:,:,0], labels=possibleParams.astype(int))
pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange,
                                               trialsEachCond[:,:,0], labels=freqLabels)
plt.setp(pRaster, ms=markerSize)
extraplots.set_ticks_fontsize(axFreqTuningNoLaser, fontSizeTicks)
axFreqTuningNoLaser.set_xticklabels([])
plt.ylabel('Laser OFF\nFreq. (kHz)', fontsize=fontSizeLabel)
#plt.title('Freq tuning')

# -------------------------------------------------------------------------------

plt.sca(axFreqTuningLaser)
pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial, timeRange,
                                               trialsEachCond[:,:,1], labels=freqLabels)
plt.setp(pRaster, ms=markerSize, color=laserColor)
plt.xlabel('Time (s)', fontsize=fontSizeLabel)
extraplots.set_ticks_fontsize(axFreqTuningLaser, fontSizeTicks)
plt.ylabel('Laser ON\nFreq. (kHz)', fontsize=fontSizeLabel)

# -------------------------------------------------------------------------------

plt.sca(axFreqTuningBoth)
logPossibleParams = np.log2(possibleParams)
plt.plot(logPossibleParams, meanFiringRateEachCond[:,0],'o-', color='k', mec='none')
plt.plot(logPossibleParams, meanFiringRateEachCond[:,1],'o-', color=laserColor, mec='none')
#axFreqTuningBoth.set_xscale('log')
#xtickLabels = [str(int(2**x)) for x in axFreqTuningBoth.get_xticks()]
xTicksHz = [2000, 4000, 8000, 16000, 32000]
xTicks = np.log2(xTicksHz)
axFreqTuningBoth.set_xticks(xTicks)
xtickLabels = [str(int(x/1000)) for x in xTicksHz]
axFreqTuningBoth.set_xticklabels(xtickLabels)
extraplots.boxoff(axFreqTuningBoth)
extraplots.set_ticks_fontsize(axFreqTuningBoth, fontSizeTicks)
plt.ylabel('Firing rate (spk/s)', fontsize=fontSizeLabel)
plt.xlabel('Frequency (kHz)', fontsize=fontSizeLabel)

plt.show()


# size = 3.6, 1.75
#plt.savefig('/tmp/terminal_inactivation_ephys.png', format='png')
figSize = [7, 3]
extraplots.save_figure('terminal_inactivation_ephys', 'svg', figSize, outputDir='/tmp/')
