"""
Estimate size of spikes and background signal for cells included in the study.

Check dates of recordings:
celldb.groupby('date')['subject'].count()

It takes ~2min to run (loading raw data from SATA HD).

"""

import os
import sys
sys.path.append('..')
import numpy as np
import json
from matplotlib import pyplot as plt
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import extraplots
import studyparams
import studyutils

SAVE = 0

# -- Data paths --
processedDataDir = '/data/neuropixels'
rawDataDir = '/media/sjara/My Book/neuropixels/'
# WARNING!: sessions for acid006 on 2023-03-22 were stored under the name inpi001
#           These sessions are not currently processed by this script.

# -- Load cell database --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_freqtuning.h5')
celldbAll = celldatabase.load_hdf(dbFilename)

# -- Select subset of cells (copied from figure_overall_firing.py) --
maxChangeFactor = studyparams.MAX_CHANGE_FACTOR
responsive = studyutils.find_tone_responsive_cells(celldbAll, frThreshold=5)
steadyParams = ['ToneBaselineFiringRate', 'ToneAvgEvokedFiringRate'] 
steady = studyutils.find_steady_cells(celldbAll, steadyParams, maxChangeFactor)
includedCells = responsive & steady
celldb = celldbAll[includedCells].copy()

# -- Prepare additional columns --
celldb['amplitudeFactor'] = np.nan
celldb['bitVolts'] = np.nan
celldb['meanAmplitude'] = np.nan
celldb['stdAmplitude'] = np.nan
celldb['minAmplitude'] = np.nan
celldb['maxAmplitude'] = np.nan
celldb['templatePeakToPeak'] = np.nan
celldb['templateExtreme'] = np.nan
celldb['stDevTrace'] = np.nan
celldb['rawDataExists'] = False

# -- Load raw traces --
#celldb2 = celldb[celldb.index==636] # 617, 636
for indIter, (indRow, dbRow) in enumerate(celldb.iterrows()):
    subject = dbRow['subject']
    sessionDate = dbRow['date']
    sessionTime = dbRow['ephysTime'][0]
    pdepth = dbRow['pdepth']
    session = f'{sessionDate}_{sessionTime}'
    clusterID = dbRow['cluster']
    bestChan = dbRow['bestChannel']
    
    multiSessionDir = os.path.join(processedDataDir, subject, f'multisession_{sessionDate}_{pdepth}um_processed')
    singleSessionDir = os.path.join(processedDataDir, subject, f'{session}_processed_multi')
    recDataDir = os.path.join(rawDataDir, subject, session, 'Record Node 101/experiment1/recording1/')
    rawDataFile = os.path.join(recDataDir, 'continuous/Neuropix-PXI-100.0/continuous.dat')
    structFile = os.path.join(recDataDir, 'structure.oebin')

    if not os.path.exists(structFile) and sessionDate == '2023-03-22' and subject == 'acid006':
        print(f'{indIter} [{indRow}] File not found: {structFile}')
        #acid006/2023-03-22_13-51-20 # Does not exist
        #acid006/2023-05-11_12-45-45 # use instead
        session = '2023-05-11_12-45-45'
        recDataDir = os.path.join(rawDataDir, subject, session, 'Record Node 101/experiment1/recording1/')
        rawDataFile = os.path.join(recDataDir, 'continuous/Neuropix-PXI-100.0/continuous.dat')
        structFile = os.path.join(recDataDir, 'structure.oebin')
        celldb.at[indRow, 'rawDataExists'] = False
    else:
        celldb.at[indRow, 'rawDataExists'] = True
        
    with open(structFile, 'r') as file:
        stfile = json.load(file)
    print(f'{indIter} [{indRow}] Loaded structure file: {structFile}')
    nChannels = stfile['continuous'][0]['num_channels']
    sampleRate = stfile['continuous'][0]['sample_rate']
    bitVolts = stfile['continuous'][0]['channels'][0]['bit_volts']
    rawdata = np.memmap(rawDataFile, dtype='int16', mode='c')
    samplesPerChan = rawdata.shape[0]//nChannels
    rawdata = rawdata.reshape((samplesPerChan, nChannels))
    nSamplesToProcess = 60 * sampleRate  # 1 min
    traceToProcess = rawdata[:nSamplesToProcess, bestChan] * bitVolts
    stDevTrace = np.std(traceToProcess)
    celldb.at[indRow, 'stDevTrace'] = stDevTrace
        
    '''
    # -- Load raw data --
    if os.path.exists(structFile):
        with open(structFile, 'r') as file:
            stfile = json.load(file)
        print(f'{indIter} [{indRow}] Loaded structure file: {structFile}')
        nChannels = stfile['continuous'][0]['num_channels']
        sampleRate = stfile['continuous'][0]['sample_rate']
        bitVolts = stfile['continuous'][0]['channels'][0]['bit_volts']
        rawdata = np.memmap(rawDataFile, dtype='int16', mode='c')
        samplesPerChan = rawdata.shape[0]//nChannels
        rawdata = rawdata.reshape((samplesPerChan, nChannels))
        nSamplesToProcess = 60 * sampleRate  # 1 min
        traceToProcess = rawdata[:nSamplesToProcess, bestChan] * bitVolts
        stDevTrace = np.std(traceToProcess)
        celldb.at[indRow, 'stDevTrace'] = stDevTrace
        celldb.at[indRow, 'rawDataExists'] = True
    else:
        # Hardcoded values from acid010 2023-08-07_10-53-44
        sampleRate = 30000
        bitVolts = 0.19499999284744262695  
        print(f'{indIter} [{indRow}] File not found: {structFile}')
        #continue
    '''
    
    # -- Load processed data --
    amplitudeAll = np.load(os.path.join(multiSessionDir, 'amplitudes.npy')).squeeze()
    clusters = np.load(os.path.join(singleSessionDir, 'spike_clusters.npy')).squeeze()
    spikeTimesAll = np.load(os.path.join(singleSessionDir, 'spike_times.npy')).squeeze()
    #templates = np.load(os.path.join(multiSessionDir, 'templates.npy'))
    #(nOrigClusters, nTimePoints, nChannels) = templates.shape
    templatesBestChan = np.load(os.path.join(multiSessionDir, 'spike_shapes.npy'))
    
    spikeInds = np.flatnonzero(clusters==clusterID)
    nSpikes = len(spikeInds)
    spikeTimes = spikeTimesAll[spikeInds]
    amplitude = amplitudeAll[spikeInds]

    #templateBestChan = templates[clusterID, :, bestChan]  # When loading templates
    templateBestChan = templatesBestChan[clusterID, :]  # When loading spike_shapes
    #templateMin = np.abs(np.min(templateBestChan))
    templateExtreme = np.max(np.abs(templateBestChan))
    templateNormalized = templateBestChan/templateExtreme
    templatePeakToPeak = np.ptp(templateNormalized)
    amplitudeFactor = 180*bitVolts*templateExtreme  # Estimated from the traces
    scaledAmplitude = amplitudeFactor * amplitude
    meanAmplitude = scaledAmplitude.mean()
    stdAmplitude = scaledAmplitude.std()
    minAmplitude = scaledAmplitude.min()
    maxAmplitude = scaledAmplitude.max()

    celldb.at[indRow, 'amplitudeFactor'] = amplitudeFactor
    celldb.at[indRow, 'bitVolts'] = bitVolts
    celldb.at[indRow, 'meanAmplitude'] = meanAmplitude
    celldb.at[indRow, 'stdAmplitude'] = stdAmplitude
    celldb.at[indRow, 'minAmplitude'] = minAmplitude
    celldb.at[indRow, 'maxAmplitude'] = maxAmplitude
    celldb.at[indRow, 'templatePeakToPeak'] = templatePeakToPeak
    celldb.at[indRow, 'templateExtreme'] = templateExtreme

# -- Plot results --
plt.clf()
#plt.hist(celldb['meanAmplitude'], bins=40)
plt.hist(celldb.stDevTrace[celldb.rawDataExists], 11, fc='0.7', ec='w', alpha=1);
plt.hist(celldb['meanAmplitude'], bins=40, fc='darkolivegreen', ec='w', alpha=1)
plt.xlabel('Mean amplitude or StDev of trace (uV)')
plt.ylabel('N cells')
plt.ylim(0, 80)
plt.xlim(0, 400)
extraplots.boxoff(plt.gca())
plt.legend(['Voltage trace Std Dev', 'Mean spike amplitude'])
plt.show()
#extraplots.save_figure('temp_spikes_snr', 'svg', [6, 2], '/tmp/')

if SAVE:
    dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
    outputFilename = os.path.join(dbPath,f'celldb_{studyparams.STUDY_NAME}_spikeSNR.h5')
    celldatabase.save_hdf(celldb, outputFilename)
