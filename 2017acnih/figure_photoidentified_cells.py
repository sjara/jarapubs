''' 
Create figure showing bandwidth tuning of photoidentified cells.

Laser trials: at least 100 trials of 100ms laser, average iti 0.8 seconds
Bandwidth tuning: 30 trials each bandwidth, sound 1 sec long, average iti 1.5 seconds
Center frequency determined with shortened tuning curve (16 freq, two intensities, best frequency one that elicits most spikes)
AM rate selected as one eliciting highest spike rate and most consistent response
'''
import os
import sys
import numpy as np
from matplotlib import pyplot as plt
from jaratoolbox import colorpalette as cp
from jaratoolbox import extraplots
reload(extraplots)
from jaratoolbox import settings
import matplotlib.gridspec as gridspec
import matplotlib
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none'


FIGNAME = 'photoidentified_cells_bandwidth_tuning'
dataDir = os.path.join(settings.FIGURES_DATA_PATH, '2017acnih', FIGNAME)

PANELS_TO_PLOT = [1,1,1]  # [Laser, ExperimentalTuning, ModelTuning]

'''
# -- Define which excitatory cell to use --
excitatoryCells = ['band016_2016-12-11_T6_c6.npz','band002_2016-08-11_T4_c5.npz','band002_2016-08-12_T6_c4.npz','band003_2016-08-18_T6_c6.npz']
args = sys.argv[1:]
if len(args):
    cellToUse = excitatoryCells[int(args[0])]
else:
    cellToUse = excitatoryCells[0]
print cellToUse
#cellFileNames = ['band004_2016-09-09_T6_c4.npz', 'band015_2016-11-12_T8_c4.npz', cellToUse]
'''

# Old version: PV, SOM, Exc
#cellFileNames = ['band004_2016-08-30_T4_c5.npz', 'band015_2016-11-12_T8_c4.npz', 'band016_2016-12-11_T6_c6.npz']
#cellColor = [cp.TangoPalette['Chameleon3'], cp.TangoPalette['ScarletRed1'], cp.TangoPalette['SkyBlue2']]

#cellFileNames = ['band016_2016-12-11_T6_c6.npz', 'band015_2016-11-12_T8_c4.npz', 'band004_2016-08-30_T4_c5.npz' ]
cellFileNames = ['band002_2016-08-12_T6_c4.npz', 'band015_2016-11-12_T8_c4.npz', 'band004_2016-08-30_T4_c5.npz' ]
cellColor = [cp.TangoPalette['SkyBlue2'], cp.TangoPalette['ScarletRed1'], cp.TangoPalette['Chameleon3'] ]


SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'photoidentified_bandwidth_tuning' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [14,6]

fontSizeLabels = 14 #12
fontSizeTicks = 12 #10
fontSizePanel = 16
labelDis = 0.1
labelPosX = [0.017, 0.23, 0.61]   # Horiz position for panel labels  0.44,
labelPosY = [0.96, 0.45]    # Vert position for panel labels
laserColor = 'c'

fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')


# -- Load laser example --
gs0 = gridspec.GridSpec(4,3)
gs0.update(top=0.9, left=0.05, right=0.7, wspace=0.6, hspace=0.1)

# --- Raster plot of laser response ---
if PANELS_TO_PLOT[0]:
    indc = 2
    laserFilename = 'example_laser_response_'+cellFileNames[indc]
    laserDataFullPath = os.path.join(dataDir,laserFilename)
    laserData = np.load(laserDataFullPath)
    laserDuration = 0.1
    
    axLaser = plt.subplot(gs0[3,0])
    plt.cla()
    axLaser.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
                     fontsize=fontSizePanel, fontweight='bold')
    axLaser.annotate('B', xy=(labelPosX[0],0.32), xycoords='figure fraction',
                     fontsize=fontSizePanel, fontweight='bold')
    nTrials = 30
    laserIndexLimitsEachTrial = laserData['indexLimitsEachTrial'][:,:nTrials]
    laserTimeRange = [-0.1,0.3]
    pRaster, hcond, zline = extraplots.raster_plot(laserData['spikeTimesFromEventOnset'],
                                                   laserIndexLimitsEachTrial,
                                                   laserTimeRange)
    axLaser.set_xticks([-0.1, 0, 0.1, 0.2, 0.3])
    axLaser.set_yticks([0,nTrials])
    yLims = np.array(plt.ylim())
    laserBarBottom = 1.1*yLims[1]
    laserBarTop = laserBarBottom + 0.1*(yLims[1]-yLims[0])
    plt.hold(1)
    plt.fill([0,laserDuration,laserDuration,0],[laserBarBottom,laserBarBottom,laserBarTop,laserBarTop],
             fc=laserColor, ec='none', clip_on=False)
    #plt.plot([0,laserDuration],laserBarBottom*[1,1], color=laserColor, clip_on=False)
    plt.hold(0)
    plt.setp(pRaster, ms=3, color=cellColor[indc])
    plt.setp(hcond, visible=False)
    plt.setp(zline, visible=False)
    extraplots.boxoff(axLaser)
    extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
    plt.ylabel('Trials',fontsize=fontSizeLabels)
    plt.xlabel('Time from laser onset (s)',fontsize=fontSizeLabels)

    
    #extraplots.breakaxis(0.15, plt.gca().get_ylim()[0], np.diff(plt.gca().get_xlim())/20.0, np.diff(plt.gca().get_ylim())/10.0)




    
gs1 = gridspec.GridSpec(3,5)
gs1.update(top=0.95, left=0.1, right=0.98, wspace=0.35, hspace=0.2)


# -- Plots for each cell type --
if PANELS_TO_PLOT[1]:
    for indc,cell in enumerate(cellFileNames):

        bandFilename = 'example_bandwidth_tuning_'+cell
        bandDataFullPath = os.path.join(dataDir,bandFilename)
        bandData = np.load(bandDataFullPath)


        # --- Raster plot of sound response at different bandwidths ---
        axRaster = plt.subplot(gs1[indc,1])
        axRaster.annotate('C', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')
        timeRange = [-0.2,1.3]
        ####colorsEachCond = 
        pRaster, hcond, zline = extraplots.raster_plot(bandData['spikeTimesFromEventOnset'],
                                                       bandData['indexLimitsEachTrial'],
                                                       timeRange,
                                                       trialsEachCond=bandData['trialsEachCond'][:,:,-1],
                                                       labels=bandData['firstSortLabels'])
        plt.setp(pRaster, ms=2, color=cellColor[indc])
        axRaster.set_xticks([0,0.5,1])
        extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
        plt.ylabel('Bandwidth (oct)',fontsize=fontSizeLabels)
        if indc==2:
            axRaster.set_xlabel('Time from sound onset (s)',fontsize=fontSizeLabels)
        else:
            axRaster.set_xticklabels('')


        # --- Plot of bandwidth tuning ---
        spikeArray = bandData['spikeArray'][:,-1].flatten()
        errorArray = bandData['errorArray'][:,-1].flatten()
        bands = bandData['possibleBands']
        plt.subplot(gs1[indc,2])
        plt.plot(range(len(bands)), spikeArray, '-o', ms=7, lw=3,
                 color=cellColor[indc], mec=cellColor[indc], clip_on=False)
        plt.fill_between(range(len(bands)), spikeArray - errorArray, 
                         spikeArray + errorArray, alpha=0.2, color='0.5', edgecolor='none')
        axCurve = plt.gca()
        #axLaser.annotate('D', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction',
        #                 fontsize=fontSizePanel, fontweight='bold')
        axCurve.set_xticklabels(bands)
        plt.ylabel('Firing rate (spk/s)',fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
        if indc==2:
            axCurve.set_xlabel('Bandwidth (oct)',fontsize=fontSizeLabels)
        else:
            axCurve.set_xticklabels('')
        extraplots.boxoff(axCurve)
        xLims = plt.xlim(); yLims = plt.ylim()
        extraplots.breakaxis(0.5, yLims[0], np.diff(xLims)/40.0, np.diff(yLims)/20.0, gap=0.5)
        extraplots.breakaxis(5.5, yLims[0], np.diff(xLims)/40.0, np.diff(yLims)/20.0, gap=0.5)
        plt.xlim(xLims); plt.ylim(yLims)
        if indc==0:
            plt.title('Mouse A1',fontsize=fontSizeLabels,fontweight='normal')


# -- Plot model curves --
modelDataDir = './modeldata'
if PANELS_TO_PLOT[2] & os.path.isdir(modelDataDir):
    import pandas as pd
    #modelDataFiles = ['SSNbandwidthTuning_noRFwidth_regime1.csv','SSNbandwidthTuning_noRFwidth_regime2.csv']
    modelDataFiles = ['SSNbandwidthTuning_NewOct42017__regime1_Fig3.csv','SSNbandwidthTuning_noRFwidth_regime2.csv']
    titleStrings = ['Model 1', 'Model 2']
    for indm, oneModelFile in enumerate(modelDataFiles):
        modelData = pd.read_csv(os.path.join(modelDataDir,oneModelFile))
        modelBW = modelData['BW(oct)']
        modelRates = [modelData['y_PV'], modelData['y_SOM'], modelData['y_E']]

        for indc,rates in enumerate(modelRates):
            axModel = plt.subplot(gs1[indc,3+indm])
            #plt.plot(np.log2(modelBW), rates, 'o', lw=5, color=cellColor[indc], mec=cellColor[indc])
            #plt.plot(modelBW[1:], rates[1:], 'o-', lw=5, color=cellColor[indc], mec=cellColor[indc], clip_on=True)
            plt.plot(np.log2(modelBW[1:]), rates[1:], '-', lw=5, color=cellColor[indc], mec=cellColor[indc], clip_on=True)
            plt.xlim(np.log2([1./8, 8]))
            #xTicks = axModel.get_xticks()
            xTicks = np.log2([1./4, 1./2, 1, 2, 4])
            axModel.set_xticks(xTicks)
            newTickLabels = [str(val) for val in (2**np.array(xTicks))]
            axModel.set_xticklabels(newTickLabels)
            #axModel.set_xticklabels(bands)
            plt.ylabel('Firing rate (spk/s)',fontsize=fontSizeLabels)
            extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
            if indc==2:
                axModel.set_xlabel('Bandwidth (oct)',fontsize=fontSizeLabels)
            else:
                axModel.set_xticklabels('')
            if indc==0:
                plt.title(titleStrings[indm],fontsize=fontSizeLabels,fontweight='normal')
            extraplots.boxoff(axModel)
        axModel.annotate('D', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')

plt.show()



if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
