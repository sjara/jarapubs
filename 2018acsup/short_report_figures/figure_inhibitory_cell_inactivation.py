import os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.colors
import matplotlib.patches as patches

from scipy import stats

from jaratoolbox import settings
from jaratoolbox import extraplots


import figparams
reload(figparams)


FIGNAME = 'figure_inhibitory_cell_inactivation'
dataDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)

PANELS = [1,1,1,1,1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'Fig2_effect_of_inhibitory_inactivation_on_suppression' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figFormat = 'svg'
figSize = [9,10] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
fontSizeLegend = figparams.fontSizeLegend

labelPosX = [0.01, 0.3, 0.63, 0.37, 0.61]   # Horiz position for panel labels
labelPosY = [0.98, 0.66, 0.32]    # Vert position for panel labels

ExcColour = figparams.colp['excitatoryCell']
Exclight = matplotlib.colors.colorConverter.to_rgba(ExcColour, alpha=0.5)
PVcolour = figparams.colp['PVcell']
PVlight = matplotlib.colors.colorConverter.to_rgba(PVcolour, alpha=0.5)
SOMcolour = figparams.colp['SOMcell']
SOMlight = matplotlib.colors.colorConverter.to_rgba(SOMcolour, alpha=0.5)

laserColour = figparams.colp['greenLaser']

soundColour = figparams.colp['sound']

fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gs = gridspec.GridSpec(3,3,width_ratios=[0.8,0.8,1.1], height_ratios=[1,1,0.9])
gs.update(top=0.98, left=0.1, right=0.95, bottom=0.08, wspace=0.4, hspace=0.3)

#exampleNoSOM = 'example_SOM_inactivation_band055_2018-03-16_1000um_T2_c4.npz'
#exampleNoSOM = 'example_SOM_inactivation_band055_2018-03-16_1400um_T4_c6.npz'
exampleNoSOM = 'example_SOM_inactivation_band055_2018-03-20_1200um_T6_c4.npz'
#exampleNoSOM = 'example_SOM_inactivation_band073_2018-09-14_1200um_T1_c6.npz'

#exampleNoPV = 'example_PV_inactivation_band056_2018-03-23_1300um_T2_c4.npz'
exampleNoPV = 'example_PV_inactivation_band062_2018-05-24_1300um_T2_c2.npz'
#exampleNoPV = 'example_PV_inactivation_band062_2018-05-25_1250um_T4_c2.npz'

summaryFileName = 'all_inactivated_cells_stats.npz'


def bootstrap_median_CI(data, reps=1000, interval=95):
    medians = np.zeros(reps)
    for ind in range(reps):
        samples = np.random.choice(data, len(data), replace=True)
        medians[ind] = np.median(samples)
    low = np.percentile(medians,(100-interval)/2.0)
    high = np.percentile(medians,interval+(100-interval)/2.0)
    return [low, high]
        
    

# --- Raster plots of sound response with and without laser ---
if PANELS[0]:
    exampleCells = [exampleNoPV, exampleNoSOM]
    
    cellColours = [[PVcolour, PVlight],[SOMcolour, SOMlight]]
    
    panelLabels = ['b', 'e']
    
    for indCell, cell in enumerate(exampleCells):
        axRaster = gs[indCell, 1]
        
        exampleDataFullPath = os.path.join(dataDir,cell)
        exampleData = np.load(exampleDataFullPath)
        colours = [[ExcColour,Exclight], cellColours[indCell]]
        possibleBands = exampleData['possibleBands']
        bandLabels = possibleBands.tolist()
        bandLabels[-1] = 'WN'
        
        laserTrials = exampleData['possibleLasers']
        inner = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=axRaster, wspace=0.1, hspace=0.25)
        for laser in laserTrials:
            thisAx = plt.subplot(inner[laser,0])
            #colorEachCond = np.tile(colours[laser], len(possibleBands)/2+1)
            colorEachCond = colours[laser]*(len(possibleBands)/2+1)
            pRaster, hcond, zline = extraplots.raster_plot(exampleData['spikeTimesFromEventOnset'],
                                                       exampleData['indexLimitsEachTrial'],
                                                       exampleData['timeRange'],
                                                       trialsEachCond=exampleData['trialsEachCond'][:,1:,laser],
                                                       labels=bandLabels[1:],
                                                       colorEachCond=colorEachCond)
            plt.setp(pRaster, ms=3)
            plt.ylabel('Bandwidth (oct)', fontsize=fontSizeLabels)
            yLims = np.array(plt.ylim())
            soundPatch = patches.Rectangle((0.0,yLims[1]*1.03),1.0,yLims[1]*0.05,linewidth=1,edgecolor=soundColour,facecolor=soundColour,clip_on=False)
            thisAx.add_patch(soundPatch)
            if not laser:
                thisAx.set_xticklabels('')
            else:
                laserPatch = patches.Rectangle((-0.1,yLims[1]*1.11),1.3,yLims[1]*0.05,linewidth=1,edgecolor=laserColour,facecolor=laserColour,clip_on=False)
                thisAx.add_patch(laserPatch)
        plt.xlabel('Time from sound onset (s)', fontsize=fontSizeLabels)
        thisAx.annotate(panelLabels[indCell], xy=(labelPosX[1],labelPosY[indCell]), xycoords='figure fraction',
                             fontsize=fontSizePanel, fontweight='bold')
        
    #also do the panel labels for the cartoons in here I guess
    cartoonPanels = ['a', 'd']
    for indPanel, panel in enumerate(cartoonPanels):
        thisAx.annotate(panel, xy=(labelPosX[0],labelPosY[indPanel]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')


# --- Plot of bandwidth tuning with and without laser ---
if PANELS[1]:
    exampleCells = [exampleNoPV, exampleNoSOM]
    figLegends = ['no PV+', 'no SOM+']
    
    excColour = figparams.colp['excitatoryCell']
    PVcolour = figparams.colp['PVcell']
    SOMcolour = figparams.colp['SOMcell']
    colours = [[excColour, PVcolour],[excColour, SOMcolour]]
    lineType = ['-', '--']
        
    panelLabels = ['c', 'f']
    legendLocs = ['upper right', 'center right']
    
    for indCell, cell in enumerate(exampleCells):
        exampleDataFullPath = os.path.join(dataDir,cell)
        exampleData = np.load(exampleDataFullPath)
        
        sustainedResponseArray = exampleData['sustainedResponseArray']
        sustainedSEM = exampleData['sustainedSEM']
        baseline = sustainedResponseArray[0]
    
        bands = exampleData['possibleBands']
        laserTrials = exampleData['possibleLasers']
        bandLabels = ['{}'.format(band) for band in np.unique(possibleBands)]
        
        fitBands = exampleData['fitBandsNoZero']
        fitResponseNoLaser = exampleData['fitResponseNoZeroNoLaser']
        fitResponseLaser = exampleData['fitResponseNoZeroLaser']
        fits = [fitResponseNoLaser, fitResponseLaser]
        
        lines = []
        SIs = [float(exampleData['suppIndNoZeroNoLaser']), float(exampleData['suppIndNoZeroLaser'])]
        print "{} SIs: {}".format(figLegends[indCell], SIs)
        
        plt.hold(1)
        axCurve = plt.subplot(gs[indCell,2])
        axCurve.set_xscale('log', basex=2)
        for laser in laserTrials:
            thisResponse = sustainedResponseArray[:,laser].flatten()
            thisSEM = sustainedSEM[:,laser].flatten()
            plt.plot(bands[1:], thisResponse[1:], 'o', ms=5,
                     color=colours[indCell][laser], mec=colours[indCell][laser], clip_on=False)
            plt.errorbar(bands[1:], thisResponse[1:], yerr = [thisSEM[1:], thisSEM[1:]], 
                         fmt='none', ecolor=colours[indCell][laser])
            line, = plt.plot(fitBands, fits[laser], lineType[laser], lw=1.5, color=colours[indCell][laser])
            lines.append(line)
        #plt.legend([lines[-1],lines[0]],['{}, SI = {:.2f}'.format(figLegends[indCell], SIs[1]),'Control, SI = {:.2f}'.format(SIs[0])], loc='best', frameon=False, fontsize=fontSizeLabels)
        plt.legend([lines[-1],lines[0]],[figLegends[indCell],'control'], loc=legendLocs[indCell], frameon=False, fontsize=fontSizeLegend, borderpad=2.0)
        axCurve.annotate(panelLabels[indCell], xy=(labelPosX[2],labelPosY[indCell]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')
        axCurve.set_ylim(bottom=0)
        axCurve.set_xticks(bands[1:])
        axCurve.tick_params(top=False, right=False, which='both')
        extraplots.boxoff(axCurve)
        extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
#         if not indCell:
#             axCurve.set_xticklabels('')
#         if indCell:
            #axCurve.set_xticks(bands)
        bands = bands.tolist()
        bands[-1] = 'WN'
        #bands[1::2] = ['']*len(bands[1::2]) #remove every other label so x axis less crowded
        axCurve.set_xticklabels(bands[1:])
        #plt.xticks(rotation=-45)
        plt.xlabel('Bandwidth (oct)',fontsize=fontSizeLabels)
        plt.ylabel('Firing rate (spk/s)',fontsize=fontSizeLabels)
        plt.xlim(0.2,7) #expand x axis so you don't have dots on y axis

# bottom three panels have their own gridspec of sorts
axSummaries = gs[2,:]
gs2 = gridspec.GridSpecFromSubplotSpec(1,3, subplot_spec=axSummaries, wspace=0.5, hspace=0.2, width_ratios=[1,0.4,1])

# --- plot scatter plot of suppression with and without laser ---   
if PANELS[2]:
    summaryDataFullPath = os.path.join(dataDir,summaryFileName)
    summaryData = np.load(summaryDataFullPath)
    
    PVsustainedSuppressionNoLaser = summaryData['fitPVsustainedSuppressionNoZeroNoLaser']
    PVsustainedSuppressionLaser = summaryData['fitPVsustainedSuppressionNoZeroLaser']
    
    SOMsustainedSuppressionNoLaser = summaryData['fitSOMsustainedSuppressionNoZeroNoLaser']
    SOMsustainedSuppressionLaser = summaryData['fitSOMsustainedSuppressionNoZeroLaser']

    PVsupNoLaser = summaryData['fitMeanPVsupNoZeroNoLaser']
    PVsupLaser = summaryData['fitMeanPVsupNoZeroLaser']
    
    semPVsupNoLaser = summaryData['fitSemPVsupNoZeroNoLaser']
    semPVsupLaser = summaryData['fitSemPVsupNoZeroLaser']
    
    SOMsupNoLaser = summaryData['fitMeanSOMsupNoZeroNoLaser']
    SOMsupLaser = summaryData['fitMeanSOMsupNoZeroLaser']
    
    semSOMsupNoLaser = summaryData['fitSemSOMsupNoZeroNoLaser']
    semSOMsupLaser = summaryData['fitSemSOMsupNoZeroLaser']
    
    panelLabel = 'g'
    
    cellLabels = ['no PV+', 'no SOM+']
    cellTypeColours = [PVcolour, SOMcolour]
    
    axScatter = plt.subplot(gs2[0,0])
    
    plt.hold(True)
    plt.plot([-2,2],[-2,2], 'k--')
    l2, = plt.plot(SOMsustainedSuppressionNoLaser,SOMsustainedSuppressionLaser, 'o', color='none', mec=SOMcolour, ms=3.2, markeredgewidth=1.2)
    l1, = plt.plot(PVsustainedSuppressionNoLaser,PVsustainedSuppressionLaser, 's', color=PVcolour, mec='none', ms=4)
    plt.ylabel('Suppression Index (laser)',fontsize=fontSizeLabels)
    plt.xlabel('Suppression Index (control)',fontsize=fontSizeLabels)
    plt.xlim(-0.1,1.1)
    plt.ylim(-0.1,1.1)

    plt.legend([l1,l2], cellLabels, loc='upper left', fontsize=fontSizeLegend, numpoints=1, handlelength=0.3, markerscale=1.7, frameon=False,)

    extraplots.boxoff(axScatter)
    axScatter.set(adjustable='box-forced', aspect='equal')
    
    axScatter.annotate(panelLabel, xy=(labelPosX[0],labelPosY[2]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')
    
    # compute stats and print in terminal
    noPV = stats.wilcoxon(PVsustainedSuppressionNoLaser,PVsustainedSuppressionLaser)[1]
    noSOM = stats.wilcoxon(SOMsustainedSuppressionNoLaser, SOMsustainedSuppressionLaser)[1]
    print "Change in SI p values:\nno PV: {0}\nno SOM: {1}".format(noPV,noSOM)
    
# --- plot median change in suppression for each cell type inactivated ---
if PANELS[3]:
    summaryDataFullPath = os.path.join(dataDir,summaryFileName)
    summaryData = np.load(summaryDataFullPath)
    
    PVsustainedSuppressionNoLaser = summaryData['fitPVsustainedSuppressionNoZeroNoLaser']
    PVsustainedSuppressionLaser = summaryData['fitPVsustainedSuppressionNoZeroLaser']
    
    SOMsustainedSuppressionNoLaser = summaryData['fitSOMsustainedSuppressionNoZeroNoLaser']
    SOMsustainedSuppressionLaser = summaryData['fitSOMsustainedSuppressionNoZeroLaser']
    
    noPVsupDiff = PVsustainedSuppressionLaser-PVsustainedSuppressionNoLaser
    noSOMsupDiff = SOMsustainedSuppressionLaser-SOMsustainedSuppressionNoLaser
    
    panelLabel = 'h'
    
    cellLabels = ['no PV+', 'no SOM+']
    cellTypeColours = [PVcolour, SOMcolour]
    
    supDiffs = [noPVsupDiff, noSOMsupDiff]
    
    axBox = plt.subplot(gs2[0,1])
    
    width = 0.7
    loc = np.arange(1,3)
#     SIMeans = [np.mean(noPVsupDiff), np.mean(noSOMsupDiff)]
#     SISEMs = [stats.sem(noPVsupDiff), stats.sem(noSOMsupDiff)]
#     for indType in range(len(SIMeans)):
#         plt.plot([loc[indType]-width/2,loc[indType]+width/2], [SIMeans[indType],SIMeans[indType]], color=cellTypeColours[indType], linewidth=3) #means
#         plt.errorbar(loc[indType], SIMeans[indType], yerr = SISEMs[indType], 
#                          fmt='none', ecolor=cellTypeColours[indType], lw=1.5, capsize=5)
#         
#     plt.plot([-5,5], [0,0], 'k--', zorder=-10)

    SIMedians = [np.median(noPVsupDiff), np.median(noSOMsupDiff)]
    SICIs = [bootstrap_median_CI(noPVsupDiff),bootstrap_median_CI(noSOMsupDiff)]
    for indType in range(len(SIMedians)): 
        plt.plot([loc[indType]-width/2,loc[indType]+width/2], [SIMedians[indType],SIMedians[indType]], color=cellTypeColours[indType], linewidth=3) #medians
        # MAKING THE ERROR BARS MANUALLY BECAUSE plt.errorbars WAS TOO MUCH A PAIN IN THE ASS
        plt.plot([loc[indType], loc[indType]],SICIs[indType], color=cellTypeColours[indType], linewidth=1.5) #error bars
        plt.plot([loc[indType]-width/8,loc[indType]+width/8],[SICIs[indType][0],SICIs[indType][0]], color=cellTypeColours[indType], linewidth=1.5) #bottom caps
        plt.plot([loc[indType]-width/8,loc[indType]+width/8],[SICIs[indType][1],SICIs[indType][1]], color=cellTypeColours[indType], linewidth=1.5) #top caps
     
    plt.plot([-5,5], [0,0], 'k--', zorder=-10)
    
#     plt.hold(True)
#     bplot = plt.boxplot(supDiffs, widths=0.6, showfliers=False)
#        
#     for box in range(len(bplot['boxes'])):
#         plt.setp(bplot['boxes'][box], color=cellTypeColours[box], linewidth=2)
#         plt.setp(bplot['whiskers'][2*box:2*(box+1)], color='none')
#         plt.setp(bplot['caps'][2*box:2*(box+1)], color='none')
#         plt.setp(bplot['medians'][box], color='k', linewidth=3)
#    
#     plt.setp(bplot['medians'], color='k')
    axBox.set_xticks([1,2])
    axBox.set_xticklabels(cellLabels)

    #plt.plot([-5,5], [0,0], 'k-', zorder=-10)
    plt.xlim(0.4, 2.6)
    axBox.set_xticks(range(1,len(supDiffs)+1))
    axBox.set_xticklabels(cellLabels, fontsize=fontSizeLabels, rotation=-45)
    yLims = (-0.2,0.1)
    #yLims = (-0.15, 0.1)
    plt.ylim(yLims)
    extraplots.significance_stars([1,2], yLims[1]*0.92, yLims[1]*0.07, gapFactor=0.2)
    extraplots.boxoff(axBox)
    plt.ylabel('Change in Suppression Index',fontsize=fontSizeLabels)
    
    axBox.annotate(panelLabel, xy=(labelPosX[3],labelPosY[2]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')
    
    pVal = stats.ranksums(noPVsupDiff, noSOMsupDiff)[1]
    print "PV vs SOM change in SI p val: {}".format(pVal)

# --- plot peak vs white noise change in firing rate with inactivation ---    
if PANELS[4]:
    summaryDataFullPath = os.path.join(dataDir,summaryFileName)
    summaryData = np.load(summaryDataFullPath)
    
    PVpeakChangeFR = summaryData['fitPVpeakChangeFRNoZero']
    PVWNChangeFR = summaryData['fitPVWNChangeFRNoZero']
    
    SOMpeakChangeFR = summaryData['fitSOMpeakChangeFRNoZero']
    SOMWNChangeFR = summaryData['fitSOMWNChangeFRNoZero']
    
    PVpeakChange = summaryData['fitMeanPVpeakChangeNoZero']
    PVWNChange = summaryData['fitMeanPVWNChangeNoZero']
    
    semPVpeakChange = summaryData['fitSemPVpeakChangeNoZero']
    semPVWNChange = summaryData['fitSemPVWNChangeNoZero']
    
    SOMpeakChange = summaryData['fitMeanSOMpeakChangeNoZero']
    SOMWNChange = summaryData['fitMeanSOMWNChangeNoZero']
    
    semSOMpeakChange = summaryData['fitSemSOMpeakChangeNoZero']
    semSOMWNChange = summaryData['fitSemSOMWNChangeNoZero']
    
    panelLabel = 'i'
    
    cellLabels = ['no PV+', 'no SOM+']
    
    PVcolour = figparams.colp['PVcell']
    SOMcolour = figparams.colp['SOMcell']
    cellTypeColours = [PVcolour, SOMcolour]
    
    axScatter = plt.subplot(gs2[0,2])
    
    plt.hold(True)
    plt.plot([-20,30],[-20,30], 'k--')
    l2, = plt.plot(SOMpeakChangeFR,SOMWNChangeFR, 'o', color='none', mec=SOMcolour, ms=3.2, markeredgewidth=1.2)
    l1, = plt.plot(PVpeakChangeFR,PVWNChangeFR, 's', color=PVcolour, mec='none', ms=4)
    plt.ylabel('Change in response to WN (spk/s)',fontsize=fontSizeLabels)
    plt.xlabel('Change in response to \n preferred bandwidth (spk/s)',fontsize=fontSizeLabels)
    plt.xlim(-5,8)
    plt.ylim(-5,8)
#     plt.xlim(-10,22)
#     plt.ylim(-10,22)
    plt.legend([l1,l2], cellLabels, loc='upper left', fontsize=fontSizeLegend, numpoints=1, handlelength=0.3, markerscale=1.7, frameon=False,)
    axScatter.annotate(panelLabel, xy=(labelPosX[4],labelPosY[2]), xycoords='figure fraction',
                         fontsize=fontSizePanel, fontweight='bold')
    extraplots.boxoff(axScatter)
    axScatter.set(adjustable='box-forced', aspect='equal')
    
    # compute stats and print in terminal
    noPV = stats.wilcoxon(PVpeakChangeFR,PVWNChangeFR)[1]
    noSOM = stats.wilcoxon(SOMpeakChangeFR, SOMWNChangeFR)[1]
    print "Change in peak  vs WN response p values:\nno PV: {0}\nno SOM: {1}".format(noPV,noSOM)
    
    
if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)