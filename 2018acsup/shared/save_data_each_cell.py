"""
Save the data for each cell separately.

It takes about 4 seconds.
"""

import os
import numpy as np
from jaratoolbox import celldatabase
from jaratoolbox import ephyscore

dbPath = '/data/figuresdata/2018acsup/shared/celldb_2018acsup_photoid.h5'
outputDataDir = '/data/2018acsup_photoid'
dbPathNew = os.path.join(outputDataDir,'dataframe_2018acsup_photoid_reduced.h5')

celldb = celldatabase.load_hdf(dbPath)

FIRING_RATE_THRESHOLD = 2  # spk/s
celldb = celldb.query('baseFiringRate50ms>@FIRING_RATE_THRESHOLD or ' +
                      'laserFiringRate50ms>@FIRING_RATE_THRESHOLD')

for indCell, (indRow, dbRow) in enumerate(celldb.iterrows()):
    oneCell = ephyscore.Cell(dbRow, useModifiedClusters=True)
    ephysData, noBehav = oneCell.load('laserPulse')

    spikeTimes = ephysData['spikeTimes']
    laserOnsetTimes = ephysData['events']['laserOn']
    #laserOffsetTimes = ephysData['events']['laserOff']

    cellStr = str(oneCell).replace(' ','_')
    filename = os.path.join(outputDataDir, f'{cellStr}.npz')

    np.savez(filename, indRow=indRow, spikeTimes=spikeTimes, laserOnsetTimes=laserOnsetTimes)
    print(f'{indCell}/{len(celldb)} Saved {filename}')
    
# -- Create simpler database --
celldbLaser = celldb.copy()
###celldbLaser.rename(columns={"genotype": "mouseGenotype"}, inplace=True)

celldbLaser.drop(columns=['behavSuffix', 'info', 'paradigm', 'recordingTrack',
                          'sessionType'], inplace=True)

# -- Add columns with cell identity --
# This follows the parameters specified in studyparams.py (line 20 & 33-35)
LASER_RESPONSE_PVAL = 0.001 # We want to be EXTRA sure not to include false positives
EXC_LASER_RESPONSE_PVAL = 0.5 #for selecting putative excitatory cells NOT responsive to laser
EXC_SPIKE_WIDTH = 0.0004
spikePeakTimes = np.array(list(celldbLaser['spikePeakTimes']))
celldbLaser['spikeWidth'] = spikePeakTimes[:,2]-spikePeakTimes[:,1]

celldbLaser['identity10ms'] = 'undef'
PV10ms = ((celldbLaser.mouseGenotype == 'PV::ChR2') &
          (celldbLaser.pValLaser10ms < LASER_RESPONSE_PVAL) &
          (celldbLaser.laserFiringRate10ms > celldbLaser.baseFiringRate10ms))
SOM10ms = ((celldbLaser.mouseGenotype == 'SOM::ChR2') &
           (celldbLaser.pValLaser10ms < LASER_RESPONSE_PVAL) &
           (celldbLaser.laserFiringRate10ms > celldbLaser.baseFiringRate10ms))
EXC10ms = ((celldbLaser.mouseGenotype == 'SOM::ChR2') &
           (celldbLaser.pValLaser10ms > EXC_LASER_RESPONSE_PVAL) &
           (celldbLaser.laserFiringRate10ms <= celldbLaser.baseFiringRate10ms) &
           (celldbLaser.spikeWidth > EXC_SPIKE_WIDTH))
celldbLaser.loc[PV10ms,'identity10ms'] = 'PV'
celldbLaser.loc[SOM10ms,'identity10ms'] = 'SOM'
celldbLaser.loc[EXC10ms,'identity10ms'] = 'EXC'

celldbLaser['identity50ms'] = 'undef'
PV50ms = ((celldbLaser.mouseGenotype == 'PV::ChR2') &
          (celldbLaser.pValLaser50ms < LASER_RESPONSE_PVAL) &
          (celldbLaser.laserFiringRate50ms > celldbLaser.baseFiringRate50ms))
SOM50ms = ((celldbLaser.mouseGenotype == 'SOM::ChR2') &
           (celldbLaser.pValLaser50ms < LASER_RESPONSE_PVAL) &
           (celldbLaser.laserFiringRate50ms > celldbLaser.baseFiringRate50ms))
EXC50ms = ((celldbLaser.mouseGenotype == 'SOM::ChR2') &
           (celldbLaser.pValLaser50ms > EXC_LASER_RESPONSE_PVAL) &
           (celldbLaser.laserFiringRate50ms <= celldbLaser.baseFiringRate50ms) &
           (celldbLaser.spikeWidth > EXC_SPIKE_WIDTH))
celldbLaser.loc[PV50ms,'identity50ms'] = 'PV'
celldbLaser.loc[SOM50ms,'identity50ms'] = 'SOM'
celldbLaser.loc[EXC50ms,'identity50ms'] = 'EXC'

celldbLaser.to_hdf(dbPathNew, key='celldb', mode='w') 

