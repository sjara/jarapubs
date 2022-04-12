'''
For reward size modulation during 0.05-0.35sec from center-out.
Generate and store intermediate data for plots showing center-out-aligned firing activity of ac neurons recorded in reward-change task (for trials using the low and high frequency separately). 
For raster data, output contains spikeTimestamps, eventOnsetTimes, spikeTimesFromEventOnset, condEachTrial, labelEachTrial, as well as meta params.
For psth data, output contains spikeCountMat, timeVec, condEachTrial, as well as meta params.

Lan Guo 20171027
'''
import os
import sys
import numpy as np
import pandas as pd
from jaratoolbox import loadbehavior
from jaratoolbox import loadopenephys
from jaratoolbox import spikesanalysis
from jaratoolbox import behavioranalysis
from jaratoolbox import settings

STUDY_NAME = '2017rc'

FIGNAME = 'reward_modulation_movement_selective_cells'
dataDir = os.path.join(settings.FIGURES_DATA_PATH, STUDY_NAME, FIGNAME)

if not os.path.exists(dataDir):
    os.mkdir(dataDir)

colorDict = {'leftMoreLowFreq':'g',
             'rightMoreLowFreq':'m',
             'sameRewardLowFreq':'y',
             'leftMoreHighFreq':'r',
             'rightMoreHighFreq':'b',
             'sameRewardHighFreq':'darkgrey'}

# -- These example cells I picked manually  --#
cellParamsList = []

exampleCell = {'subject':'gosi001',
              'date':'2017-05-06',
              'tetrode':3,
               'cluster':5,
               'brainRegion':'ac'} # high freq
cellParamsList.append(exampleCell)

exampleCell = {'subject':'gosi004',
              'date':'2017-02-13',
              'tetrode':7,
               'cluster':8,
               'brainRegion':'ac'} # low freq
cellParamsList.append(exampleCell)

exampleCell = {'subject':'gosi004',
              'date':'2017-03-15',
              'tetrode':6,
               'cluster':7,
               'brainRegion':'ac'} # high freq
cellParamsList.append(exampleCell)

exampleCell = {'subject':'gosi004',
              'date':'2017-03-18',
              'tetrode':3,
               'cluster':4,
               'brainRegion':'ac'} # high freq
cellParamsList.append(exampleCell)

exampleCell = {'subject':'gosi004',
              'date':'2017-04-06',
              'tetrode':7,
               'cluster':5,
               'brainRegion':'ac'} # high freq
cellParamsList.append(exampleCell)

exampleCell = {'subject':'gosi008',
              'date':'2017-03-14',
              'tetrode':7,
               'cluster':8,
               'brainRegion':'ac'} # high freq
cellParamsList.append(exampleCell)

# -- Here we can choose to generate data for a specific cell instead of every cell -- #
cellIndToGenerate = False

####################################################################################
scriptFullPath = os.path.realpath(__file__)
timeRange = [-0.5,1]
binWidth = 0.010
EPHYS_SAMPLING_RATE = 30000.0
soundTriggerChannel = 0
paradigm = '2afc'
minBlockSize = 30
freqsToPlot = ['low', 'high']
alignment = 'center-out'
###################################################################################
def get_trials_each_cond_reward_change(bdata, freqToPlot, colorCondDict, byBlock=True, minBlockSize=30):
    '''Function to generate selection vector showing which trials to plot for each behavior conditions and the color to use in the plot label.
    :param arg1: bdata object (with missing trials removed).
    :param arg2: A string indicating which frequency to plot, value of 'low' or 'high'.
    :param arg3: Boolean indicating whether to split the plot by behavior blocks. 
    :param arg4: An integer for the minimal number of trials a block has to have to be plotted.
    :param arg5: A dictionary indicating which color label each behavior condition gets in the final plot.
    
    :return: trialsEachCond, colorsEachCond, labelEachCond
    '''
    # -- Select trials to plot from behavior file -- #
    correct = bdata['outcome']==bdata.labels['outcome']['correct']
    possibleFreq = np.unique(bdata['targetFrequency'])
    numFreqs = len(possibleFreq)
    currentBlock = bdata['currentBlock']
    # -- Select trials to plot based on desired frequencies to plot and whether to plot by block -- #
    ### Recordings during reward change usually have 2 frequencies, low freq means go to left, right freq means go to right ###
    if freqToPlot == 'low':
        freq = possibleFreq[0] 

    elif freqToPlot == 'high':
        freq = possibleFreq[-1]

    oneFreq = bdata['targetFrequency'] == freq #vector for selecing trials presenting this frequency
    correctOneFreq = oneFreq  & correct 

    # -- Find trials each block (if plotting by block) or find trials each type (e.g. more_left, more_right; if not plotting by block) -- #
    if byBlock:
        bdata.find_trials_each_block()
        numBlocks = bdata.blocks['nBlocks']
        trialsEachBlock = bdata.blocks['trialsEachBlock']
        correctTrialsEachBlock = trialsEachBlock & correctOneFreq[:,np.newaxis]
        correctBlockSizes = sum(correctTrialsEachBlock)
        if (correctBlockSizes[-1] < minBlockSize): #A check to see if last block is too small to plot
            correctTrialsEachBlock = correctTrialsEachBlock[:,:-1]
            numBlocks -= 1
        trialsEachCond = correctTrialsEachBlock
        
        colorEachCond = np.empty(numBlocks, dtype=object)
        labelEachCond = np.empty(numBlocks, dtype=object)
        #pdb.set_trace()
        for blockNum in range(numBlocks):
            currentBlockLabel = currentBlock[trialsEachBlock[:,blockNum]][0]
            if freqToPlot == 'low':
                if currentBlockLabel == bdata.labels['currentBlock']['same_reward']:
                    colorEachCond[blockNum] = colorCondDict['sameRewardLowFreq']
                    labelEachCond[blockNum] = 'low freq same reward'
                elif currentBlockLabel == bdata.labels['currentBlock']['more_left']:
                    colorEachCond[blockNum] = colorCondDict['leftMoreLowFreq'] 
                    labelEachCond[blockNum] = 'low freq left more'
                elif currentBlockLabel == bdata.labels['currentBlock']['more_right']:   
                    colorEachCond[blockNum] = colorCondDict['rightMoreLowFreq']
                    labelEachCond[blockNum] = 'low freq right more'
            elif freqToPlot == 'high':
                if currentBlockLabel == bdata.labels['currentBlock']['same_reward']:
                    colorEachCond[blockNum] = colorCondDict['sameRewardHighFreq']
                    labelEachCond[blockNum] = 'high freq same reward'
                elif currentBlockLabel == bdata.labels['currentBlock']['more_left']:
                    colorEachCond[blockNum] = colorCondDict['leftMoreHighFreq'] 
                    labelEachCond[blockNum] = 'high freq left more'
                elif currentBlockLabel == bdata.labels['currentBlock']['more_right']:   
                    colorEachCond[blockNum] = colorCondDict['rightMoreHighFreq']
                    labelEachCond[blockNum] = 'high freq right more'

    else:
        blockTypes = [bdata.labels['currentBlock']['same_reward'],bdata.labels['currentBlock']['more_left'],bdata.labels['currentBlock']['more_right']]
        trialsEachType = behavioranalysis.find_trials_each_type(currentBlock,blockTypes)
        oneFreqCorrectBlockSameReward = correctOneFreq&trialsEachType[:,0]
        oneFreqCorrectBlockMoreLeft = correctOneFreq&trialsEachType[:,1]
        oneFreqCorrectBlockMoreRight = correctOneFreq&trialsEachType[:,2]
        
        trialsEachCond = np.c_[oneFreqCorrectBlockSameReward,oneFreqCorrectBlockMoreLeft,oneFreqCorrectBlockMoreRight]
        if freqToPlot == 'low':
            colorEachCond = [colorCondDict['sameRewardLowFreq'],colorCondDict['leftMoreLowFreq'],colorCondDict['rightMoreLowFreq']]
            labelEachCond = ['low freq same reward', 'low freq left more', 'low freq right more']
        elif freqToPlot == 'high':
            colorEachCond = [colorCondDict['sameRewardHighFreq'],colorCondDict['leftMoreHighFreq'],colorCondDict['rightMoreHighFreq']]
            labelEachCond = ['high freq same reward', 'high freq left more', 'high freq right more']
    return trialsEachCond, colorEachCond, labelEachCond

###################################################################################
# -- Access mounted behavior and ephys drives for psycurve and switching mice -- #
BEHAVIOR_PATH = settings.BEHAVIOR_PATH_REMOTE
EPHYS_PATH = settings.EPHYS_PATH_REMOTE

if not os.path.ismount(BEHAVIOR_PATH):
    os.system('sshfs -o idmap=user jarauser@jarahub:/data/behavior/ {}'.format(BEHAVIOR_PATH))

if not os.path.ismount(EPHYS_PATH):
    os.system('sshfs -o idmap=user jarauser@jarastore:/data2016/ephys/ {}'.format(EPHYS_PATH))


# -- Select an example cell from allcells file -- #
if cellIndToGenerate:
    cellParamsList = [cellParamsList[cellIndToGenerate]]
for cellParams in cellParamsList:
    animal = cellParams['subject']
    date = cellParams['date']
    tetrode = cellParams['tetrode']
    cluster = cellParams['cluster']
    brainRegion = cellParams['brainRegion']
    celldbPath = os.path.join(settings.DATABASE_PATH, '{}_database.h5'.format(animal))
    celldb = pd.read_hdf(celldbPath, key='reward_change')
    
    ### Using cellDB methode to find this cell in the cellDB ###
    oneCell = celldb.loc[(celldb.subject==animal) & (celldb.date==date) & (celldb.tetrode==tetrode) & (celldb.cluster==cluster)]
    sessionsThisCell = oneCell.iloc[0].sessiontype
    rcInd = sessionsThisCell.index('behavior')
    rcEphysThisCell = oneCell['ephys'].iloc[0][rcInd]
    rcBehavThisCell = oneCell['behavior'].iloc[0][rcInd]

    ## Get behavior data associated with 2afc session ###
    behavFileName = rcBehavThisCell
    behavFile = os.path.join(BEHAVIOR_PATH,animal,behavFileName)
    bdata = loadbehavior.FlexCategBehaviorData(behavFile,readmode='full')


    ### Get events data ###
    fullEventFilename=os.path.join(EPHYS_PATH, animal, rcEphysThisCell, 'all_channels.events')
    eventData = loadopenephys.Events(fullEventFilename)
    ##### Get event onset times #####
    eventData.timestamps = np.array(eventData.timestamps)/EPHYS_SAMPLING_RATE #hard-coded ephys sampling rate!!


    ### GEt spike data of just this cluster ###
    spikeFilename = os.path.join(EPHYS_PATH, animal, rcEphysThisCell, 'Tetrode{}.spikes'.format(tetrode))
    spikeData = loadopenephys.DataSpikes(spikeFilename)
    spikeData.timestamps = spikeData.timestamps/EPHYS_SAMPLING_RATE
    clustersDir = os.path.join(EPHYS_PATH, animal, rcEphysThisCell)+'_kk'
    clusterFilename = os.path.join(clustersDir, 'Tetrode{}.clu.1'.format(tetrode))
    clusters = np.fromfile(clusterFilename, dtype='int32', sep=' ')[1:]
    spikeData.timestamps = spikeData.timestamps[clusters==cluster]
    spikeData.samples = spikeData.samples[clusters==cluster, :, :]
    spikeData.samples = spikeData.samples.astype(float)-2**15# FIXME: this is specific to OpenEphys
    # FIXME: This assumes the gain is the same for all channels and records
    spikeData.samples = (1000.0/spikeData.gain[0,0]) * spikeData.samples
    #spikeData = ephyscore.CellData(oneCell) #This defaults to settings ephys path
    spikeTimestamps = spikeData.timestamps

    # -- Check to see if ephys has skipped trials, if so remove trials from behav data -- #
    eventOnsetTimes=np.array(eventData.timestamps)
    soundOnsetEvents = (eventData.eventID==1) & (eventData.eventChannel==soundTriggerChannel)
    soundOnsetTimeEphys = eventOnsetTimes[soundOnsetEvents]
    soundOnsetTimeBehav = bdata['timeTarget']

    # Find missing trials
    missingTrials = behavioranalysis.find_missing_trials(soundOnsetTimeEphys,soundOnsetTimeBehav)
    # Remove missing trials
    bdata.remove_trials(missingTrials)

    # -- calculate response and reaction time, as well as trials each type -- #
    responseTimes = bdata['timeSideIn'] - bdata['timeCenterOut']
    reactionTimes = bdata['timeCenterOut'] - bdata['timeTarget']
    currentBlock = bdata['currentBlock']
    blockTypeLabels = ['more_left', 'more_right']
    blockTypes = [bdata.labels['currentBlock']['more_left'],bdata.labels['currentBlock']['more_right']]
    trialsEachType = behavioranalysis.find_trials_each_type(currentBlock,blockTypes)
    possibleFreq = np.unique(bdata['targetFrequency'])
    correct = bdata['outcome']==bdata.labels['outcome']['correct']

    # -- Select trials to plot from behavior file -- #
    for freq in freqsToPlot:
        trialsEachCond, colorEachCond, labelEachCond = get_trials_each_cond_reward_change(bdata, freqToPlot=freq, byBlock=True, minBlockSize=30, colorCondDict=colorDict)
        
        # -- Calculate eventOnsetTimes aligned to center exit -- #
        soundOnsetEvents = (eventData.eventID==1) & (eventData.eventChannel==soundTriggerChannel)
        EventOnsetTimes = eventOnsetTimes[soundOnsetEvents]
        diffTimes=bdata['timeCenterOut']-bdata['timeTarget']
        EventOnsetTimes+=diffTimes

        # -- Calculate arrays for plotting raster -- #
        (spikeTimesFromEventOnset,trialIndexForEachSpike,indexLimitsEachTrial) = \
        spikesanalysis.eventlocked_spiketimes(spikeTimestamps,EventOnsetTimes,timeRange)

        # -- Save behavior related times and labels -- #
        if freq == 'low':
            oneFreqTrials = bdata['targetFrequency'] == possibleFreq[0] 
        elif freq == 'high':
            oneFreqTrials = bdata['targetFrequency'] == possibleFreq[-1]
        leftMoreTrialsOneFreq = trialsEachType[:,0] & oneFreqTrials & correct
        rightMoreTrialsOneFreq = trialsEachType[:,1] & oneFreqTrials & correct
        outputFile = 'behavior_times_{}freq_{}_{}.npz'.format(freq, animal,date)
        outputFullPath = os.path.join(dataDir,outputFile)
        np.savez(outputFullPath, leftMoreTrials=leftMoreTrialsOneFreq, rightMoreTrials=rightMoreTrialsOneFreq, responseTimes=responseTimes, reactionTimes=reactionTimes)

        # -- Save raster intermediate data -- #    
        #outputDir = os.path.join(settings.FIGURESDATA, figparams.STUDY_NAME)
        outputFile = 'example_rc_{}aligned_raster_{}freq_{}_{}_T{}_c{}.npz'.format(alignment, freq, animal, date, tetrode, cluster)
        outputFullPath = os.path.join(dataDir,outputFile)
        np.savez(outputFullPath, spikeTimestamps=spikeTimestamps, eventOnsetTimes=EventOnsetTimes, spikeTimesFromEventOnset=spikeTimesFromEventOnset, indexLimitsEachTrial=indexLimitsEachTrial, condLabels=labelEachCond, trialsEachCond=trialsEachCond, colorEachCond=colorEachCond, script=scriptFullPath, EPHYS_SAMPLING_RATE=EPHYS_SAMPLING_RATE, soundTriggerChannel=soundTriggerChannel, timeRange=timeRange, frequencyPloted=freq, alignedTo=alignment, **cellParams)

        # -- Calculate additional arrays for plotting psth -- #
        timeVec = np.arange(timeRange[0],timeRange[-1],binWidth)
        spikeCountMat = spikesanalysis.spiketimes_to_spikecounts(spikeTimesFromEventOnset,indexLimitsEachTrial,timeVec)

        # -- Save psth intermediate data -- #
        #outputDir = os.path.join(settings.FIGURESDATA, figparams.STUDY_NAME)
        outputFile = 'example_rc_{}aligned_psth_{}freq_{}_{}_T{}_c{}.npz'.format(alignment, freq, animal, date, tetrode, cluster)
        outputFullPath = os.path.join(dataDir,outputFile)
        print 'Saving {0} ...'.format(outputFullPath)
        np.savez(outputFullPath, spikeCountMat=spikeCountMat, timeVec=timeVec, condLabels=labelEachCond, trialsEachCond=trialsEachCond, colorEachCond=colorEachCond, timeRange=timeRange, binWidth=binWidth, EPHYS_SAMPLING_RATE=EPHYS_SAMPLING_RATE, soundTriggerChannel=soundTriggerChannel, script=scriptFullPath, frequencyPloted=freq, alignedTo=alignment, **cellParams)
