"""
Show example rasters and spike shapes given cell database and files with spike times.
"""

import os
import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

dataDir = '/data/2018acsup_photoid/'
dfPath = os.path.join(dataDir,'dataframe_2018acsup_photoid_reduced.h5')
celldf = pd.read_hdf(dfPath)

# -- Create an array with the spike waveform for each cell --
spikeWaveforms = np.array(list(celldf.spikeShape))  # [nCells, nSamples]


def eventlocked_spiketimes(timeStamps, eventOnsetTimes, timeRange):
    """
    Align spike timestamps to an event.

    Args:
        timeStamps: (np.array) the time of each spike.
        eventOnsetTimes: (np.array) the time of each instance of the event to lock to.
        timeRange: (list) two-element list specifying time-range to extract around event.

    Returns:
        spikeTimesFromEventOnset: 1D array with time of spikes locked to event.
        trialIndexForEachSpike: 1D array with the trial corresponding to each spike.
           The first spike index is 0.
        indexLimitsEachTrial: [2,nTrials] range of spikes for each trial. Note that
           the range is from firstSpike to lastSpike+1 (like in python slices)
    """
    nTrials = len(eventOnsetTimes)
    spikeTimesFromEventOnset = np.empty(0,dtype='float64')
    spikeIndices = np.empty(0,dtype='int')
    trialIndexForEachSpike = np.empty(0,dtype='int')
    indexLimitsEachTrial = np.empty((2,nTrials),dtype='int')
    accumIndexFirstSpike = 0
    for indtrial in np.arange(nTrials):
        thisTrialRange = eventOnsetTimes[indtrial] + timeRange
        firstSpikeInTrial = np.searchsorted(timeStamps,thisTrialRange[0])
        # NOTE: lastSpikeInTrial must not be negative, because slice(0,-1) means from 0 to last.
        lastSpikeInTrialPlusOne = np.searchsorted(timeStamps,thisTrialRange[-1])
        spikesThisTrial = np.arange(firstSpikeInTrial,lastSpikeInTrialPlusOne)
        nSpikesThisTrial = lastSpikeInTrialPlusOne - firstSpikeInTrial

        spikeTimesFromEventOnset = np.concatenate((spikeTimesFromEventOnset,
                                        timeStamps[spikesThisTrial]-eventOnsetTimes[indtrial]))
        trialIndexForEachSpike = np.concatenate((trialIndexForEachSpike,
                                            np.repeat(indtrial,nSpikesThisTrial)))
        indexLimitsEachTrial[:,indtrial] = [accumIndexFirstSpike,accumIndexFirstSpike+nSpikesThisTrial]
    return (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial)


# -- Load a few cells and plot a spike raster aligned to laser onset --
nCellsToLoad = 5
gs = gridspec.GridSpec(nCellsToLoad, 2, width_ratios=[0.3, 0.7])

timeRange = [-0.3, 0.4]  # In seconds

plt.clf()
for indCell, (indRow, dfRow) in enumerate(celldf.iterrows()):

    dataFile = (f'{dfRow.subject}_{dfRow.date}_{dfRow.pdepth}um_' +
                f'g{dfRow.egroup}c{dfRow.cluster}.npz')

    # -- Load the spikes and laser timestamps data --
    tsData = np.load(os.path.join(dataDir,dataFile))
    spikeTimes = tsData['spikeTimes']
    eventOnsetTimes = tsData['laserOnsetTimes']
    
    # -- Plot the spike waveform --
    axWaveform = plt.subplot(gs[indCell%nCellsToLoad, 0])
    assert indRow==tsData['indRow']  # Make sure the waveform and timestamps are from the same cell
    axWaveform.plot(spikeWaveforms[indCell,:], lw=4)
    axWaveform.set_ylabel(dfRow['identity50ms'], rotation=0, fontweight='bold')
    
    # -- Plot the spike raster --
    (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
        eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange)

    axRaster = plt.subplot(gs[indCell%nCellsToLoad, 1])
    axRaster.plot(spikeTimesFromEventOnset, trialIndexForEachSpike, '.k', ms=2)
    axRaster.set_xlim(timeRange)
    axRaster.set_title(dataFile)
    axRaster.set_ylabel('Trials')
    if indCell%nCellsToLoad == (nCellsToLoad-1):
        axRaster.set_xlabel('Time (s)')
        #plt.show(); plt.pause(0.5); plt.clf()
        break
    else:
        axRaster.set_xticklabels([])

plt.show()


