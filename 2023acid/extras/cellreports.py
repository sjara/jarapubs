"""
Make a page report for each cell.
"""

import os
import sys
sys.path.append('..')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import ephyscore
from jaratoolbox import spikesanalysis
from jaratoolbox import extraplots
from jaratoolbox import behavioranalysis
import studyparams
from importlib import reload
reload(studyparams)


DEBUG = 0
SAVE = 0

timeRange = {'Full':     [-0.3, 0.45],
             'Baseline': [-0.2, 0],
             'Evoked':   [0.015, 0.115] }

dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_oddball.h5')
celldb = celldatabase.load_hdf(dbFilename)

measurements = ['BaselineFiringRate', 'EvokedFiringRate']
stimConditions = ['oddball', 'standard']

if DEBUG:
    indRow = 46
    indRow = 1
    celldbToUse = celldb.iloc[[indRow]]
else:
    celldbToUse = celldb

for indRow, dbRow in celldbToUse.iterrows():
    fig = plt.gcf()
    fig.clf()
    fig.set_facecolor('w')
    fig.set_size_inches([14, 8])
    gsMain = gridspec.GridSpec(3, 4)
    gsMain.update(left=0.05, right=0.98, top=0.9, bottom=0.08, wspace=0.25, hspace=0.3)
    
    if indRow%20==0:
        print(f'{indRow}/{len(celldb)} cells analyzed')
    oneCell = ephyscore.Cell(dbRow)
    plt.suptitle(f'{oneCell} [{indRow}]', fontweight='bold')
    for indr, reagent in enumerate(studyparams.REAGENTS):
        for inds, (sessionType, sessionInfo) in enumerate(studyparams.ODDBALL_SESSION_TYPES.items()):

            # ========= THIS CODE SHOULD BE IDENTICAL TO database_oddball_analysis.py ======
            ephysData, bdata = oneCell.load(reagent+sessionType)
            spikeTimes = ephysData['spikeTimes']
            eventOnsetTimes = ephysData['events']['stimOn']

            stimEachTrial = bdata['currentStartFreq']
            nTrials = len(stimEachTrial)

            # If the ephys data is 1 more than the bdata, delete the last ephys trial.
            if len(stimEachTrial) == len(eventOnsetTimes)-1:
                eventOnsetTimes = eventOnsetTimes[:len(stimEachTrial)]
            assert len(stimEachTrial) == len(eventOnsetTimes), \
                "Number of trials in behavior and ephys do not match for {oneCell}"
            
            (spikeTimesFromEventOnset, trialIndexForEachSpike, indexLimitsEachTrial) = \
                spikesanalysis.eventlocked_spiketimes(spikeTimes, eventOnsetTimes, timeRange['Full'])

            possibleStim = np.unique(stimEachTrial)
            trialsEachCond = behavioranalysis.find_trials_each_type(stimEachTrial, possibleStim)

            # -- Get oddball and standard trials (checking things match ODDBALL_SESSION_TYPES) --
            indOddball = np.argmin(trialsEachCond.sum(axis=0))
            indEachCond = {'oddball': indOddball, 'standard': 1-indOddball}
            # -- According to currentStartFreq, this is how the stim are ordered --
            stimOrder = {'FM': ['FM_up', 'FM_down'],
                         'Chord': ['low_freq', 'high_freq']}
            # -- Sanity check --
            oddballStim = bdata.labels['oddballStim'][bdata['oddballStim'][-1]]
            stimType = bdata.labels['stimType'][bdata['stimType'][-1]]
            assert oddballStim == stimOrder[stimType][indOddball], \
                f"Oddball stim should be {stimOrder[stimType][indOddball]}"
            # ===========================================================================
            
            thisAx = plt.subplot(gsMain[indr, inds])
            pRaster, hcond, zline = extraplots.raster_plot(spikeTimesFromEventOnset, indexLimitsEachTrial,
                                                           timeRange['Full'], trialsEachCond)
            plt.setp(pRaster, ms=0.5)
            keyOdd = reagent + sessionInfo['oddball'] + 'OddballEvokedFiringRate'
            keyStd = reagent + sessionInfo['standard'] + 'StandardEvokedFiringRate'
            #keyOdd = reagent + sessionInfo['oddball'] + 'OddballBaselineFiringRate'
            #keyStd = reagent + sessionInfo['standard'] + 'StandardBaselineFiringRate'
            pOdd = dbRow[reagent + sessionInfo['oddball'] + 'OddballPval']
            pStd = dbRow[reagent + sessionInfo['standard'] + 'StandardPval']
            frOdd = f'Odd:{dbRow[keyOdd]:0.1f}'
            frStd = f'Std:{dbRow[keyStd]:0.1f}'
            plt.title(f'{reagent} {sessionType} ({frOdd}, {frStd})\np={pOdd:0.3f}, {pStd:0.3f}')
            if indr==2:
                plt.xlabel('Time from sound onset (s)')
            else:
                thisAx.set_xticklabels([])
            #sys.exit()
    plt.show()
    plt.pause(0.01);
    plt.pause(1);

    figFormat = 'png'
    figFilename = f'cell{indRow:04}_' + str(oneCell).replace(' ','_') + '.' + figFormat
    figPath = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'cell_reports',
                           figFilename)
    if SAVE:
        plt.savefig(figPath, format=figFormat)
        os.system(f'convert {figPath} PNG8:{figPath}') # Changed in indexed to reduce size
        print(f'Saved: {figPath}')
