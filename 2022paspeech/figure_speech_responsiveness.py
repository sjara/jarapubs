"""
This creates figure 2 for 2022paspeech:
 A. Cartoon of headfixed, awake mouse ephys
 B. Diagram of sound matrix
 C. Histology image of recording track
 D. Scatter plot of recording location of each cell. Add AC areas image in inkscape.
"""

import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import colorpalette as cp
import figparams
import studyparams
import studyutils
from scipy import stats
from importlib import reload

reload(figparams)
reload(studyutils)


FIGNAME = 'selectivityIndices'
SAVE_FIGURE = 1
STATSUMMARY = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'plots_speech_responsiveness' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [5.25, 6.75] # In inches

PANELS = [3,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
fontSizeTitles = figparams.fontSizeTitles
colorSpeechResp = cp.TangoPalette['ScarletRed2']
colorNotAud = cp.TangoPalette['Aluminium3']
colorSoundResp = cp.TangoPalette['SkyBlue2']
colorStimPresented = cp.TangoPalette['Butter2']


labelPosX = [0.04, 0.32, 0.68] # Horiz position for panel labels
labelPosY = [0.97, 0.77, 0.365]    # Vert position for panel labels

plt.clf()
gsMain = gridspec.GridSpec(2, 1, height_ratios=[0.15, 0.85])
gsMain.update(left=0.07, right=0.99, top=0.98, bottom=0.07, wspace=0.35, hspace=0.25)
gsTop = gsMain[0].subgridspec(1,3, width_ratios = [0.35, 0.3, 0.35])

axCartoon = plt.subplot(gsTop[0,0])
axCartoon.set_axis_off()

axHist = plt.subplot(gsTop[0,1])
axHist.set_axis_off()

axSoundMatrix = plt.subplot(gsTop[0,2])
#axSoundMatrix.set_axis_off()

gsBottom = gsMain[1].subgridspec(2, 2, width_ratios = [0.65, 0.35], wspace=0.1, hspace=0.25)
gsDonutsAreas = gsBottom[0,1].subgridspec(2, 2)
gsDonutsRegions = gsBottom[1,1].subgridspec(2, 2)

axCellLocs = plt.subplot(gsBottom[:,0])

axDonutAudD = plt.subplot(gsDonutsAreas[0,0])
axDonutAudP = plt.subplot(gsDonutsAreas[0,1])
axDonutAudV = plt.subplot(gsDonutsAreas[1,0])
axDonutTeA = plt.subplot(gsDonutsAreas[1,1])

axDonutDP = plt.subplot(gsDonutsRegions[0,0])
axDonutDA = plt.subplot(gsDonutsRegions[0,1])
axDonutVP = plt.subplot(gsDonutsRegions[1,0])
axDonutVA = plt.subplot(gsDonutsRegions[1,1])

# -- Move Donuts a little lower --
yoffset = -0.02
for oneax in [axDonutAudD, axDonutAudP, axDonutAudV, axDonutTeA]:
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0, axpos.y0+yoffset, axpos.width, axpos.height])
yoffset = -0.03
for oneax in [axDonutDP, axDonutDA, axDonutVP, axDonutVA]:
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0, axpos.y0+yoffset, axpos.width, axpos.height])

'''
gsBottom = gsMain[1].subgridspec(4,3, width_ratios = [0.6, 0.2, 0.2])

axCellLocs = plt.subplot(gsBottom[:,0])

axDonutAudP = plt.subplot(gsBottom[0,2])
axDonutAudD = plt.subplot(gsBottom[0,1])
axDonutAudV = plt.subplot(gsBottom[1,1])
axDonutTeA = plt.subplot(gsBottom[1,2])

axDonutDP = plt.subplot(gsBottom[2,1])
axDonutDA = plt.subplot(gsBottom[2,2])
axDonutVP = plt.subplot(gsBottom[3,1])
axDonutVA = plt.subplot(gsBottom[3,2])
'''

axCartoon.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axSoundMatrix.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axHist.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axCellLocs.annotate('D', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axDonutAudP.annotate('E', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axDonutDP.annotate('F', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


# -- Plot a grid with cell borders and labels --
plt.sca(axSoundMatrix)
corner_labels = [ ['/da/', '/ta/'], ['/ba/', '/pa/'] ]
VOTstep = 20
VOTvals = np.arange(0, 60+VOTstep, VOTstep)
FTstep = 6
FTvals = np.arange(-9, 9+FTstep, FTstep)
for indFT, thisFT in enumerate(FTvals):
    for indVOT, thisVOT in enumerate(VOTvals):
        if indFT in [0, 3] or indVOT in [0, 3]:
            facecolor = colorStimPresented
        else:
            facecolor = 'none'
        onePatch = axSoundMatrix.add_patch(plt.Rectangle((thisVOT-VOTstep/2, thisFT-FTstep/2),
                                                            VOTstep, FTstep, lw=0.75, ec='black', fc=facecolor))
        onePatch.set_clip_on(False)
        if indFT in [0, 3] and indVOT in [0, 3]:
            axSoundMatrix.text(thisVOT, thisFT, corner_labels[indFT // 3][indVOT // 3],
                                ha='center', va='center', fontsize=fontSizeTicks-2)
axSoundMatrix.set_xlim(-10, 70)
axSoundMatrix.set_ylim(-12, 12)
axSoundMatrix.set_xlabel('VOT (ms)', fontsize=fontSizeLabels-2)
axSoundMatrix.set_ylabel('FT (oct/s)', fontsize=fontSizeLabels-2)
axSoundMatrix.set_xticks(VOTvals)
axSoundMatrix.set_yticks(FTvals)
axSoundMatrix.set_aspect(20/6)
extraplots.set_ticks_fontsize(axSoundMatrix, fontSizeTicks-2)
axSoundMatrix.invert_yaxis()


# -- Load responsiveness data --
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
x_coords_jittered = figData['x_coords_jittered']
z_coords_jittered = figData['z_coords_jittered']
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
isCortical = figData['isCortical']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
quantilesDV = figData['quantilesDV']
quantilesAP = figData['quantilesAP']
layersDeep = figData['layersDeep']
layer4 = figData['layer4']
layersSuperficial = figData['layersSuperficial']

'''
parser= argparse.ArgumentParser()
parser.add_argument('-d','--layerDeep', help = "only includes cells from layers 5 and 6", required = False, action='store_true')
parser.add_argument('-f','--layer4', help = "only includes cells from layer 4", required = False, action='store_true')
parser.add_argument('-s','--layerSuperficial', help = "only includes cells from layers 1 and 2/3", required = False, action='store_true')
args = parser.parse_args()
if args.layerDeep:
    useLayers = figData['layersDeep']
    print('only including layer 5 and layer 6 cells')
elif args.layer4:
    useLayers = figData['layer4']
    print('only including layer 4 cells')
elif args.layerSuperficial:
    useLayers = figData['layersSuperficial']
    print('only including layer 1 and layer 2/3 cells')
else:
    useLayers = np.ones(len(x_coords), dtype = bool)
'''
parser= argparse.ArgumentParser()
parser.add_argument('-d','--layerDeep', help = "only includes cells from layers 5 and 6", required = False, action='store_true')
parser.add_argument('-f','--layer4', help = "only includes cells from layer 4", required = False, action='store_true')
parser.add_argument('-s','--layerSuperficial', help = "only includes cells from layers 1 and 2/3", required = False, action='store_true')
args = parser.parse_args()
if args.layerDeep:
    useLayers = figData['layersDeep']
    print('only including cells in layers 5 & 6')
elif args.layer4:
    useLayers = figData['layer4']
    print('only including layer 4 cells')
elif args.layerSuperficial:
    useLayers = figData['layersSuperficial'] | figData['layer4']
    print('only including cells in layers 1-4')
else:
    useLayers = np.ones(len(x_coords), dtype = bool)


### ARE THESE USED? ###
quadrantBoundsDV = figData['quadrantBoundsDV']
quadrantBoundsAP = figData['quadrantBoundsAP']
quadrantBoundsDV_AApix = figData['quadrantBoundsDV_AApix']
quadrantBoundsAP_AApix = figData['quadrantBoundsAP_AApix']

quadrantLabels = figData['quadrantLabels']
quadrantTotals = figData['quadrantTotals']
quadrantsSpeechResponsive = figData['quadrantsSpeechResponsive']
quadrantsSoundResponsive = figData['quadrantsSoundResponsive']

boundDataFile = 'brain_areas_boundaries.npz'
boundDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
boundDataFullPath = os.path.join(boundDataDir,boundDataFile)
boundData = np.load(boundDataFullPath, allow_pickle = True)
contours = boundData['contours']
extentAP = boundData['extentAP']
extentDV = boundData['extentDV']

APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(studyutils.pix2mmAP(APtickLocs),1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round(studyutils.pix2mmDV(DVtickLocs),1)

# -- Plot the locations of the recorded cells -- #
plt.sca(axCellLocs)
studyutils.plot_quadrants(axCellLocs, extentAP, extentDV, color='0.8')
for contour in contours:
    plt.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color='k', clip_on=False)
labelSize = fontSizePanel
plt.text(234, 119, 'D', ha='center', fontsize=labelSize)
plt.text(224, 127, 'P', ha='center', fontsize=labelSize)
plt.text(233, 142, 'V', ha='center', fontsize=labelSize)
plt.text(233, 165, 'TeA', ha='center', fontsize=labelSize)

nonResp = plt.scatter(z_coords_jittered[~soundResponsive & isCortical & useLayers], y_coords[~soundResponsive & isCortical & useLayers], c = colorNotAud, s=6)
soundResp = plt.scatter(z_coords_jittered[soundResponsive & ~speechResponsive & isCortical & useLayers], y_coords[soundResponsive & ~speechResponsive & isCortical & useLayers], c = colorSoundResp, s = 6)
speechResp = plt.scatter(z_coords_jittered[speechResponsive & isCortical & useLayers] , y_coords[speechResponsive & isCortical & useLayers], c = colorSpeechResp, s=6)
plt.xlim(146, 246)
plt.ylim(220,40)
SHOW_UNITS_IN_MM = 1
if SHOW_UNITS_IN_MM:
    plt.xticks(APtickLocs, APtickLabels)
    plt.yticks(DVtickLocs, DVtickLabels)
    units = 'mm'
else:
    units = 'atlas voxels'
plt.xlabel(f'Posterior-Anterior ({units})', fontsize = fontSizeLabels)
plt.ylabel(f'Ventral-Dorsal ({units})', fontsize = fontSizeLabels)
#plt.legend([nonResp, soundResp, speechResp], ['not sound responsive', 'sound responsive','speech responsive'], loc = 'upper left', markerscale = 2 , bbox_to_anchor = (0, 1.15))
plt.legend([nonResp, soundResp, speechResp], ['Not sound responsive', 'Non-speech sound responsive','Speech responsive'], loc = 'upper center', markerscale=2, handletextpad=0.25, bbox_to_anchor = (0.555, 1.04))
axCellLocs.spines["right"].set_visible(False)
axCellLocs.spines["top"].set_visible(False)
axCellLocs.set_aspect('equal')

#plt.show(); sys.exit()


##-- Donut Plots Atlas-defined Areas--
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0]) & useLayers)
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1]) & useLayers)
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2]) & useLayers)
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3]) & useLayers)

nSoundResponsiveAudP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[0]) & useLayers)
nSoundResponsiveAudD = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[1]) & useLayers)
nSoundResponsiveAudV = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[2]) & useLayers)
nSoundResponsiveTeA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[3]) & useLayers)

nCellsAudP = np.sum((recordingAreaName == audCtxAreas[0]) & useLayers)
nCellsAudD = np.sum((recordingAreaName == audCtxAreas[1]) & useLayers)
nCellsAudV = np.sum((recordingAreaName == audCtxAreas[2]) & useLayers)
nCellsTeA = np.sum((recordingAreaName == audCtxAreas[3]) & useLayers)

nCellsDeepAudP = np.sum((recordingAreaName == audCtxAreas[0]) & layersDeep)
nCellsDeepAudD = np.sum((recordingAreaName == audCtxAreas[1]) & layersDeep)
nCellsDeepAudV = np.sum((recordingAreaName == audCtxAreas[2]) & layersDeep)
nCellsDeepTeA = np.sum((recordingAreaName == audCtxAreas[3]) & layersDeep)

nCellsSuperficialAudP = np.sum((recordingAreaName == audCtxAreas[0]) & (layersSuperficial|layer4))
nCellsSuperficialAudD = np.sum((recordingAreaName == audCtxAreas[1]) & (layersSuperficial|layer4))
nCellsSuperficialAudV = np.sum((recordingAreaName == audCtxAreas[2]) & (layersSuperficial|layer4))
nCellsSuperficialTeA = np.sum((recordingAreaName == audCtxAreas[3]) & (layersSuperficial|layer4))

titleY = 0.95
titlePad = -0

plt.sca(axDonutAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudP, nSoundResponsiveAudP - nSpeechResponsiveAudP, nCellsAudP - nSoundResponsiveAudP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudP.add_artist(circle1)
plt.title(f'AudP\n n = {nCellsAudP}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)

plt.sca(axDonutAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudD, nSoundResponsiveAudD - nSpeechResponsiveAudD, nCellsAudD - nSoundResponsiveAudD], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudD.add_artist(circle2)
plt.title(f'AudD\n n = {nCellsAudD}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)


plt.sca(axDonutAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudV, nSoundResponsiveAudV - nSpeechResponsiveAudV, nCellsAudV - nSoundResponsiveAudV], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudV.add_artist(circle3)
plt.title(f'AudV\n n = {nCellsAudV}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)

plt.sca(axDonutTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveTeA, nSoundResponsiveTeA - nSpeechResponsiveTeA, nCellsTeA - nSoundResponsiveTeA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutTeA.add_artist(circle4)
plt.title(f'TeA\n n = {nCellsTeA}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)


##-- Donut Plots Quartiles--
quadrantsUseLayers = []
quadrantsLayersDeep = []
#quadrantsLayer4 = []
quadrantsLayersSuperficial = []

a = -1
for indBinDV, thisQuantileDV in enumerate(quantilesDV):
    for indBinAP, thisQuantileAP in enumerate(quantilesAP):
        a = a+1
        quadrantsUseLayers.append(useLayers[thisQuantileDV & thisQuantileAP])
        quadrantsLayersDeep.append(layersDeep[thisQuantileDV & thisQuantileAP])
        #quadrantsLayer4.append(layer4[thisQuantileDV & thisQuantileAP])
        quadrantsLayersSuperficial.append(layersSuperficial[thisQuantileDV & thisQuantileAP]|layer4[thisQuantileDV & thisQuantileAP])

nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0] & quadrantsUseLayers[0])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1] & quadrantsUseLayers[1])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2] & quadrantsUseLayers[2])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3] & quadrantsUseLayers[3])

nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0] & quadrantsUseLayers[0])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1] & quadrantsUseLayers[1])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2] & quadrantsUseLayers[2])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3] & quadrantsUseLayers[3])

nCellsDP = quadrantTotals[0] - np.sum(~quadrantsUseLayers[0])
nCellsDA = quadrantTotals[1] - np.sum(~quadrantsUseLayers[1])
nCellsVP = quadrantTotals[2] - np.sum(~quadrantsUseLayers[2])
nCellsVA = quadrantTotals[3] - np.sum(~quadrantsUseLayers[3])




plt.sca(axDonutDP)
circle5 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveDP, nSoundResponsiveDP - nSpeechResponsiveDP, nCellsDP - nSoundResponsiveDP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutDP.add_artist(circle5)
plt.title(f'DP\n n = {nCellsDP}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)

plt.sca(axDonutDA)
circle6 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveDA, nSoundResponsiveDA - nSpeechResponsiveDA, nCellsDA - nSoundResponsiveDA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutDA.add_artist(circle6)
plt.title(f'DA\n n = {nCellsDA}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)

plt.sca(axDonutVP)
circle7 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveVP, nSoundResponsiveVP - nSpeechResponsiveVP, nCellsVP - nSoundResponsiveVP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutVP.add_artist(circle7)
plt.title(f'VP\n n = {nCellsVP}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)

plt.sca(axDonutVA)
circle8 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveVA, nSoundResponsiveVA - nSpeechResponsiveVA, nCellsVA - nSoundResponsiveVA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutVA.add_artist(circle8)
plt.title(f'VA\n n = {nCellsVA}', y=titleY, pad=titlePad, fontsize = fontSizeTicks)


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)


if STATSUMMARY:
    # Test proportion of speech responsive cells between areas (of sound responsive cells)
    oddsratio, pvalFracResponsive_AudPvsAudD = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudD],[nSoundResponsiveAudP - nSpeechResponsiveAudP, nSoundResponsiveAudD - nSpeechResponsiveAudD]]))
    oddsratio, pvalFracResponsive_AudPvsAudV = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudV],[nSoundResponsiveAudP - nSpeechResponsiveAudP, nSoundResponsiveAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudDvsAudV = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveAudV],[nSoundResponsiveAudD - nSpeechResponsiveAudD, nSoundResponsiveAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudPvsTeA = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveTeA],[nSoundResponsiveAudP - nSpeechResponsiveAudP, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudDvsTeA = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveTeA],[nSoundResponsiveAudD - nSpeechResponsiveAudD, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudVvsTeA = stats.fisher_exact(np.array([[nSpeechResponsiveAudV, nSpeechResponsiveTeA],[nSoundResponsiveAudV - nSpeechResponsiveAudV, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))

    # Test proportion of speech responsive cells between areas (of all cells)
    oddsratio, pvalFracResponsive_AudPvsAudD_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudD],[nCellsAudP - nSpeechResponsiveAudP, nCellsAudD - nSpeechResponsiveAudD]]))
    oddsratio, pvalFracResponsive_AudPvsAudV_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudV],[nCellsAudP - nSpeechResponsiveAudP, nCellsAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudDvsAudV_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveAudV],[nCellsAudD - nSpeechResponsiveAudD, nCellsAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudPvsTeA_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveTeA],[nCellsAudP - nSpeechResponsiveAudP, nCellsTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudDvsTeA_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveTeA],[nCellsAudD - nSpeechResponsiveAudD, nCellsTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudVvsTeA_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudV, nSpeechResponsiveTeA],[nCellsAudV - nSpeechResponsiveAudV, nCellsTeA - nSpeechResponsiveTeA]]))

    nQuadCompar = 6
    quadAlpha = 0.05/6
    quadComparLabels = np.array(['DP vs DA', 'DP vs. VP', 'DP vs. VA', 'DA vs. VP', 'DA vs. VA', 'VP vs. VA'])

    ##--Test frac sound responsive
    quadrantComparFracSoundResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsSoundResponsive):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracSoundResponsive = stats.fisher_exact(np.array([[np.sum(quadrantsSoundResponsive[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1])],[np.sum(quadrantsSoundResponsive[indBin]==0), np.sum(quadrantsSoundResponsive[x + indBin + 1]==0)]]))
            quadrantComparFracSoundResponsive[a] = pvalFracSoundResponsive

    ##--Test frac speech responsive out of total
    quadrantComparFracSpeechResponsive_allcells = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsSpeechResponsive):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracSpeechResponsive_allcells = stats.fisher_exact(np.array([[np.sum(quadrantsSpeechResponsive[indBin] & quadrantsUseLayers[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1] & quadrantsUseLayers[x + indBin + 1])],[len(quadrantsSoundResponsive[indBin]) -np.sum(quadrantsSpeechResponsive[indBin] & quadrantsUseLayers[indBin]), len(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsSpeechResponsive[x + indBin + 1] & quadrantsUseLayers[x + indBin + 1])]]))
            quadrantComparFracSpeechResponsive_allcells[a] = pvalFracSpeechResponsive_allcells

    ##--Test frac speech responsive out of soundResponsive
    quadrantComparFracSpeechResponsive_soundResp = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsSpeechResponsive):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracSpeechResponsive_soundResp = stats.fisher_exact(np.array([[np.sum(quadrantsSpeechResponsive[indBin]  & quadrantsUseLayers[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1] & quadrantsUseLayers[x + indBin + 1])],[np.sum(quadrantsSoundResponsive[indBin] & quadrantsUseLayers[indBin]) -np.sum(quadrantsSpeechResponsive[indBin] & quadrantsUseLayers[indBin]  & quadrantsUseLayers[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1] & quadrantsUseLayers[x + indBin + 1]) - np.sum(quadrantsSpeechResponsive[x + indBin + 1]  & quadrantsUseLayers[x + indBin + 1])]]))
            quadrantComparFracSpeechResponsive_soundResp[a] = pvalFracSpeechResponsive_soundResp

    reload(studyutils)
    # -- Print LaTeX tables --
    # NOTE: The following lines using eval() would not be necessary if the data was already
    #       stored in numpy arrays (or pandas dataframes) when initially calculated.
    brainAreas = ['AudP', 'AudD', 'AudV', 'TeA']
    nCellsEachArea = [eval('nCells'+brainArea) for brainArea in brainAreas]
    nSoundResponsiveEachArea = [eval('nSoundResponsive'+brainArea) for brainArea in brainAreas]
    nSpeechResponsiveEachArea = [eval('nSpeechResponsive'+brainArea) for brainArea in brainAreas]
    pvalsFracResponsiveAllCellsEachComp = np.full((len(brainAreas),len(brainAreas)), np.nan)
    for inda1, area1 in enumerate(brainAreas):
        for inda2, area2 in enumerate(brainAreas):
            if inda2 > inda1:
                pvalsFracResponsiveAllCellsEachComp[inda1,inda2] = eval('pvalFracResponsive_'+
                                                                       area1+'vs'+area2+'_allcells')
    pvalsFracResponsiveEachComp = np.full((len(brainAreas),len(brainAreas)), np.nan)
    for inda1, area1 in enumerate(brainAreas):
        for inda2, area2 in enumerate(brainAreas):
            if inda2 > inda1:
                pvalsFracResponsiveEachComp[inda1,inda2] = eval('pvalFracResponsive_'+area1+'vs'+area2)

    brainRegions = ['DP', 'DA', 'VP', 'VA']
    nSpeechResponsiveEachRegion = [eval('nSpeechResponsive'+brainRegion) for brainRegion in brainRegions]
    nSoundResponsiveEachRegion = [eval('nSoundResponsive'+brainRegion) for brainRegion in brainRegions]

    pvalsQuadComparFracSpeechResponsive_allcells = np.full((len(brainRegions),len(brainRegions)), np.nan)
    a = -1
    for inda1, area1 in enumerate(brainRegions):
        for inda2, area2 in enumerate(brainRegions):
            if inda2 > inda1:
                a = a+1
                pvalsQuadComparFracSpeechResponsive_allcells[inda1,inda2] = quadrantComparFracSpeechResponsive_allcells[a]


    pvalsQuadComparFracSpeechResponsive_soundResp = np.full((len(brainAreas),len(brainAreas)), np.nan)
    a = -1
    for inda1, area1 in enumerate(brainRegions):
        for inda2, area2 in enumerate(brainRegions):
            if inda2 > inda1:
                a = a+1
                pvalsQuadComparFracSpeechResponsive_soundResp[inda1,inda2] = quadrantComparFracSpeechResponsive_soundResp[a]


    table1a = studyutils.latex_table_responsive(brainAreas, nSpeechResponsiveEachArea, nCellsEachArea)
    print(table1a)
    table1b = studyutils.latex_table_pvals(brainAreas, pvalsFracResponsiveAllCellsEachComp)
    print(table1b)
    table2a = studyutils.latex_table_responsive(brainAreas,nSpeechResponsiveEachArea, nSoundResponsiveEachArea)
    print(table2a)
    table2b = studyutils.latex_table_pvals(brainAreas, pvalsFracResponsiveEachComp)
    print(table2b)
    table3a = studyutils.latex_table_responsive(brainRegions, nSpeechResponsiveEachRegion, quadrantTotals)
    print(table3a)
    table3b = studyutils.latex_table_pvals(brainRegions, pvalsQuadComparFracSpeechResponsive_allcells)
    print(table3b)
    table4a = studyutils.latex_table_responsive(brainRegions, nSpeechResponsiveEachRegion, nSoundResponsiveEachRegion)
    print(table4a)
    table4b = studyutils.latex_table_pvals(brainRegions, pvalsQuadComparFracSpeechResponsive_soundResp)
    print(table4b)



    if 1:
        # -- Print all stats --
        print('--Stats Summary--')
        print('--atlas-defined areas')
        print(f'AudP n: {nCellsAudP}, n any sound responsive: {nSoundResponsiveAudP} ({np.round((nSoundResponsiveAudP/nCellsAudP)*100, 1)}%), n speech responsive: {nSpeechResponsiveAudP}, ({np.round((nSpeechResponsiveAudP/nCellsAudP)*100, 1)}% total, {np.round((nSpeechResponsiveAudP/nSoundResponsiveAudP)*100, 1)}% soundResponsive)')
        print(f'AudD n: {nCellsAudD}, n any sound responsive: {nSoundResponsiveAudD} ({np.round((nSoundResponsiveAudD/nCellsAudD)*100, 1)}%), n speech responsive: {nSpeechResponsiveAudD}, ({np.round((nSpeechResponsiveAudD/nCellsAudD)*100, 1)}% total, {np.round((nSpeechResponsiveAudD/nSoundResponsiveAudD)*100, 1)}% soundResponsive)')
        print(f'AudV n: {nCellsAudV}, n any sound responsive: {nSoundResponsiveAudV} ({np.round((nSoundResponsiveAudV/nCellsAudV)*100, 1)}%), n speech responsive: {nSpeechResponsiveAudV}, ({np.round((nSpeechResponsiveAudV/nCellsAudV)*100, 1)}% total, {np.round((nSpeechResponsiveAudV/nSoundResponsiveAudV)*100, 1)}% soundResponsive)')
        print(f'TeA n: {nCellsTeA}, n any sound responsive: {nSoundResponsiveTeA} ({np.round((nSoundResponsiveTeA/nCellsTeA)*100, 1)}%), n speech responsive: {nSpeechResponsiveTeA}, ({np.round((nSpeechResponsiveTeA/nCellsTeA)*100, 1)}% total, {np.round((nSpeechResponsiveTeA/nSoundResponsiveTeA)*100, 1)}% soundResponsive)')
        print(f'Bonferroni corrected alpha = 0.05/6 = {np.round(0.05/6,3)}')
        print('--Frac speech responsive of sound responsive --')
        print(f'AudP vs AudD p = {np.round(pvalFracResponsive_AudPvsAudD,3)}')
        print(f'AudP vs AudV p = {np.round(pvalFracResponsive_AudPvsAudV,3)}')
        print(f'AudD vs AudV p = {np.round(pvalFracResponsive_AudDvsAudV,3)}')
        print(f'AudP vs TeA p = {np.round(pvalFracResponsive_AudPvsTeA,3)}')
        print(f'AudD vs TeA p = {np.round(pvalFracResponsive_AudDvsTeA,3)}')
        print(f'AudV vs TeA p = {np.round(pvalFracResponsive_AudVvsTeA,3)}')
        print('--Frac speech responsive of total cells --')
        print(f'AudP vs AudD p = {np.round(pvalFracResponsive_AudPvsAudD_allcells,3)}')
        print(f'AudP vs AudV p = {np.round(pvalFracResponsive_AudPvsAudV_allcells,3)}')
        print(f'AudD vs AudV p = {np.round(pvalFracResponsive_AudDvsAudV_allcells,3)}')
        print(f'AudP vs TeA p = {np.round(pvalFracResponsive_AudPvsTeA_allcells,3)}')
        print(f'AudD vs TeA p = {np.round(pvalFracResponsive_AudDvsTeA_allcells,3)}')
        print(f'AudV vs TeA p = {np.round(pvalFracResponsive_AudVvsTeA_allcells,3)}')
        print('--quadrant-defined areas')
        print(f'DP Quadrant: total n = {nCellsDP}, n soundResponsive = {nSoundResponsiveDP} ({np.round((nSoundResponsiveDP/nCellsDP)*100,1)}%), n speech responsive = {nSpeechResponsiveDP} ({np.round((nSpeechResponsiveDP/nSoundResponsiveDP)*100,1)}% of sound responsive, {np.round((nSpeechResponsiveDP/nCellsDP)*100,1)}% of total)')
        print(f'DA Quadrant: total n = {nCellsDA}, n soundResponsive = {nSoundResponsiveDA} ({np.round((nSoundResponsiveDA/nCellsDA)*100,1)}%), n speech responsive = {nSpeechResponsiveDA} ({np.round((nSpeechResponsiveDA/nSoundResponsiveDA)*100,1)}% of sound responsive, {np.round((nSpeechResponsiveDA/nCellsDA)*100,1)}% of total)')
        print(f'VP Quadrant: total n = {nCellsVP}, n soundResponsive = {nSoundResponsiveVP} ({np.round((nSoundResponsiveVP/nCellsVP)*100,1)}%), n speech responsive = {nSpeechResponsiveVP} ({np.round((nSpeechResponsiveVP/nSoundResponsiveVP)*100,1)}% of sound responsive, {np.round((nSpeechResponsiveVP/nCellsVP)*100,1)}% of total)')
        print(f'VA Quadrant: total n = {nCellsVA}, n soundResponsive = {nSoundResponsiveVA} ({np.round((nSoundResponsiveVA/nCellsVA)*100,1)}%), n speech responsive = {nSpeechResponsiveVA} ({np.round((nSpeechResponsiveVA/nSoundResponsiveVA)*100,1)}% of sound responsive, {np.round((nSpeechResponsiveVA/nCellsVA)*100,1)}% of total)')
        print('--Frac speech responsive of sound responsive --')
        print(f'DP vs DA p = {np.round(quadrantComparFracSpeechResponsive_soundResp[0],3)}')
        print(f'DP vs VP p = {np.round(quadrantComparFracSpeechResponsive_soundResp[1],3)}')
        print(f'DP vs VA p = {np.round(quadrantComparFracSpeechResponsive_soundResp[2],3)}')
        print(f'DA vs VP p = {np.round(quadrantComparFracSpeechResponsive_soundResp[3],3)}')
        print(f'DA vs VA p = {np.round(quadrantComparFracSpeechResponsive_soundResp[4],3)}')
        print(f'VP vs VA p = {np.round(quadrantComparFracSpeechResponsive_soundResp[5],3)}')
        print('--Frac speech responsive of total cells --')
        print(f'DP vs DA p = {np.round(quadrantComparFracSpeechResponsive_allcells[0],3)}')
        print(f'DP vs VP p = {np.round(quadrantComparFracSpeechResponsive_allcells[1],3)}')
        print(f'DP vs VA p = {np.round(quadrantComparFracSpeechResponsive_allcells[2],3)}')
        print(f'DA vs VP p = {np.round(quadrantComparFracSpeechResponsive_allcells[3],3)}')
        print(f'DA vs VA p = {np.round(quadrantComparFracSpeechResponsive_allcells[4],3)}')
        print(f'VP vs VA p = {np.round(quadrantComparFracSpeechResponsive_allcells[5],3)}')
        print('-- nCells in deep vs. superficial layers --')
        print(f'DP: nTotal:{len(quadrantsLayersDeep[0])}, nL1-4:{np.sum(quadrantsLayersSuperficial[0])}, nL5-6:{np.sum(quadrantsLayersDeep[0])}')
        print(f'DA: nTotal:{len(quadrantsLayersDeep[1])}, nL1-4:{np.sum(quadrantsLayersSuperficial[1])}, nL5-6:{np.sum(quadrantsLayersDeep[1])}')
        print(f'VP: nTotal:{len(quadrantsLayersDeep[2])}, nL1-4:{np.sum(quadrantsLayersSuperficial[2])}, nL5-6:{np.sum(quadrantsLayersDeep[2])}')
        print(f'VA: nTotal:{len(quadrantsLayersDeep[3])}, nL1-4:{np.sum(quadrantsLayersSuperficial[3])}, nL5-6:{np.sum(quadrantsLayersDeep[3])}')
        print(f'AudP: nTotal:{np.sum(recordingAreaName == audCtxAreas[0])}, nL1-4:{nCellsSuperficialAudP}, nL5-6:{nCellsDeepAudP}')
        print(f'AudD: nTotal:{np.sum(recordingAreaName == audCtxAreas[1])}, nL1-4:{nCellsSuperficialAudD}, nL5-6:{nCellsDeepAudD}')
        print(f'AudV: nTotal:{np.sum(recordingAreaName == audCtxAreas[2])}, nL1-4:{nCellsSuperficialAudV}, nL5-6:{nCellsDeepAudV}')
        print(f'TeA: nTotal:{np.sum(recordingAreaName == audCtxAreas[3])}, nL1-4:{nCellsSuperficialTeA}, nL5-6:{nCellsDeepTeA}')
