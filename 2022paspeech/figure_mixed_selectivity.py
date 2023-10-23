"""
This creates figure 5 (about mixed selectivity) for 2022paspeech.
 A. plot of VOT selectivity indices in A-P/D-V space. Add AC areas image in inkscape
 B. Donut plots of fraction of VOT selective cells among speech responsive cells for each Aud Area
 C. VOT selectivity indices binned in D-V space
 D. VOT selectivity indices binned in A-P space
 E-H: as in A-D, for FT
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import argparse
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import celldatabase
from jaratoolbox import colorpalette as cp
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy import stats
import studyparams
import studyutils
import figparams
from importlib import reload
reload(figparams)

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
shuffledDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
SAVE_FIGURE = 1
STATSUMMARY = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'figure_mixed_selectivity' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7.5, 5.0] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
fontSizePanel = figparams.fontSizePanel
colorMap = cm.get_cmap('Greens')
newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
colorMixedSelective = cp.TangoPalette['Orange2']
colorSingleSelective = cp.TangoPalette['SkyBlue2']
colorNotSelective = cp.TangoPalette['Aluminium3']  # 'Aluminium2'
colorBounds = cp.TangoPalette['Aluminium3'] #'0.75'


figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
z_coords_jittered = figData['z_coords_jittered']
x_coords_jittered = figData['x_coords_jittered']
subjects = np.unique(figData['subject'])
speechResponsive = figData['speechResponsive']
isCortical = figData['isCortical']
excludeSpeech = figData['excludeSpeech']
bestSelectivityIndexVot = figData['bestSelectivityIndexVot']
bestSelectivityIndexFt = figData['bestSelectivityIndexFt']
pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
shuffledVotBest = figData['shuffledVotBest']
shuffledFtBest = figData['shuffledFtBest']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
avgShuffledSIVot = np.mean(shuffledVotBest[~excludeSpeech & speechResponsive],1)
avgShuffledSIFt = np.mean(shuffledFtBest[~excludeSpeech & speechResponsive], 1)
quantilesDV = figData['quantilesDV']
quantilesAP = figData['quantilesAP']
quadrantBoundsDV = figData['quadrantBoundsDV']
quadrantBoundsAP = figData['quadrantBoundsAP']
#quadrantBoundsDV_AApix = figData['quadrantBoundsDV_AApix']
#quadrantBoundsAP_AApix = figData['quadrantBoundsAP_AApix']
quadrantLabels = figData['quadrantLabels']
quadrantTotals = figData['quadrantTotals']
quadrantsVOT = figData['quadrantsVOT']
quadrantsFT = figData['quadrantsFT']
quadrantsVotSelective = figData['quadrantsVotSelective']
quadrantsFtSelective = figData['quadrantsFtSelective']
quadrantsMixedSelective = figData['quadrantsMixedSelective']
quadrantsSingleSelective = figData['quadrantsSingleSelective']
quadrantsSpeechResponsive = figData['quadrantsSpeechResponsive']
quadrantsSoundResponsive = figData['quadrantsSoundResponsive']
quadrantTotalsByAnimal = figData['quadrantTotalsByAnimal']
quadrantsSoundResponsiveByAnimal = figData['quadrantsSoundResponsiveByAnimal']
quadrantsSpeechResponsiveByAnimal = figData['quadrantsSpeechResponsiveByAnimal']
quadrantsVotSelectiveByAnimal = figData['quadrantsVotSelectiveByAnimal']
quadrantsFtSelectiveByAnimal = figData['quadrantsFtSelectiveByAnimal']
quadrantsSingleSelectiveByAnimal = figData['quadrantsSingleSelectiveByAnimal']
quadrantsMixedSelectiveByAnimal = figData['quadrantsMixedSelectiveByAnimal']

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



boundDataFile = 'brain_areas_boundaries.npz'
boundDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
boundDataFullPath = os.path.join(boundDataDir,boundDataFile)
boundData = np.load(boundDataFullPath, allow_pickle = True)
contours = boundData['contours']
extentAP = boundData['extentAP']
extentDV = boundData['extentDV']

plt.clf()
gsMain = gridspec.GridSpec(4, 3, width_ratios = [0.40, 0.3, 0.3])
gsMain.update(left=0.07, right=0.98, top=0.97, bottom=0.05, wspace=0.3, hspace=0.4)
axMixeSelMap = plt.subplot(gsMain[:,0])
axQuadSummary = gsMain[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45], hspace=0.3)
axQuadPcts = plt.subplot(axQuadSummary[0,0])
axByAnimal = plt.subplot(axQuadSummary[1,0])
axMixDonuts = gsMain[0:,1].subgridspec(4, 1, hspace = 0.4)
axMixAudP = plt.subplot(axMixDonuts[1,0])
axMixAudD = plt.subplot(axMixDonuts[0,0])
axMixAudV = plt.subplot(axMixDonuts[2,0])
axMixTeA = plt.subplot(axMixDonuts[3,0])
# -- Move Donuts a little to the left and closer together --
xoffset = -0.03
yoffset = 0.02
for indax, oneax in enumerate([axMixAudD, axMixAudP, axMixAudV, axMixTeA]):
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0+xoffset, axpos.y0+(indax-2)*yoffset, axpos.width, axpos.height])
# -- Move third column a little higher --
yoffset = 0.04
for indax, oneax in enumerate([axQuadPcts, axByAnimal]):
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0, axpos.y0+yoffset, axpos.width, axpos.height])

#plt.subplots_adjust(top = 0.9, bottom = 0.1, hspace = 0.45, left = 0.05)

votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive & useLayers
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive  & useLayers
singleSelective = np.logical_xor(votSelective, ftSelective)
mixedSelective = votSelective & ftSelective


APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(studyutils.pix2mmAP(APtickLocs),1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round(studyutils.pix2mmDV(DVtickLocs),1)
'''
APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
'''

plt.sca(axMixeSelMap)
respNonSel = plt.scatter(z_coords_jittered[speechResponsive & ~mixedSelective & ~singleSelective & isCortical  & useLayers], y_coords[speechResponsive & ~mixedSelective & ~singleSelective & isCortical  & useLayers], c = colorNotSelective, s = 6)
singSel = plt.scatter(z_coords_jittered[speechResponsive & singleSelective & isCortical], y_coords[speechResponsive & singleSelective & isCortical], c = colorSingleSelective, s = 6)
mixSel = plt.scatter(z_coords_jittered[speechResponsive & mixedSelective & isCortical], y_coords[speechResponsive & mixedSelective & isCortical], c = colorMixedSelective, s = 6)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)
plt.legend(handles = [singSel, mixSel, respNonSel], loc="upper left",
           labels=['Single-selective', 'Mixed-selective', 'Speech-responsive, non-selective'],
           markerscale=3, handletextpad=0.25, fontsize=fontSizeTicks, bbox_to_anchor=(0.01, 1.07))
#plt.title('Mixed selectivity', fontsize = fontSizeTitles)
axMixeSelMap.set_aspect('equal')
axMixeSelMap.spines["right"].set_visible(False)
axMixeSelMap.spines["top"].set_visible(False)
axMixeSelMap.set_facecolor('none')
studyutils.plot_quadrants(axMixeSelMap, extentAP, extentDV, color=colorBounds)
for contour in contours:
    plt.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color=colorBounds, clip_on=False, zorder=-1)
labelSize = fontSizePanel
plt.text(234, 119, 'D', ha='center', fontsize=labelSize, color=colorBounds)
plt.text(224, 127, 'P', ha='center', fontsize=labelSize, color=colorBounds)
plt.text(233, 142, 'V', ha='center', fontsize=labelSize, color=colorBounds)
plt.text(233, 165, 'TeA', ha='center', fontsize=labelSize, color=colorBounds)


plt.sca(axQuadPcts)
quadrantsUseLayers = []
quadrantsMixedSelective = []
quadrantsSingleSelective = []

a = -1
for indBinDV, thisQuantileDV in enumerate(quantilesDV):
    for indBinAP, thisQuantileAP in enumerate(quantilesAP):
        a = a+1
        quadrantsUseLayers.append(useLayers[thisQuantileDV & thisQuantileAP])
        quadrantsMixedSelective.append(mixedSelective[thisQuantileDV & thisQuantileAP])
        quadrantsSingleSelective.append(singleSelective[thisQuantileDV & thisQuantileAP])



nMixedSelectiveDP = np.sum(quadrantsMixedSelective[0])
nSpeechSelectiveDP = np.sum(quadrantsMixedSelective[0]) + np.sum(quadrantsSingleSelective[0])
nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0] & quadrantsUseLayers[0])
nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0]& quadrantsUseLayers[0])
nTotalDP = quadrantTotals[0] - np.sum(~quadrantsUseLayers[0])
#fracMixedSelectiveDP = nMixedSelectiveDP/nSpeechSelectiveDP
fracMixedSelectiveDP = nMixedSelectiveDP/nSpeechResponsiveDP
DPalphaMix = 1

nMixedSelectiveDA = np.sum(quadrantsMixedSelective[1])
nSpeechSelectiveDA = np.sum(quadrantsMixedSelective[1]) + np.sum(quadrantsSingleSelective[1])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1] & quadrantsUseLayers[1])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1]& quadrantsUseLayers[1])
nTotalDA = quadrantTotals[1] - np.sum(~quadrantsUseLayers[1])
#fracMixedSelectiveDA = nMixedSelectiveDA/nSpeechSelectiveDA
fracMixedSelectiveDA = nMixedSelectiveDA/nSpeechResponsiveDA
DAalphaMix = fracMixedSelectiveDA/fracMixedSelectiveDP

nMixedSelectiveVP = np.sum(quadrantsMixedSelective[2])
nSpeechSelectiveVP = np.sum(quadrantsMixedSelective[2]) + np.sum(quadrantsSingleSelective[2])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2] & quadrantsUseLayers[2])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2] & quadrantsUseLayers[2])
nTotalVP = quadrantTotals[2] - np.sum(~quadrantsUseLayers[2])
#fracMixedSelectiveVP = nMixedSelectiveVP/nSpeechSelectiveVP
fracMixedSelectiveVP = nMixedSelectiveVP/nSpeechResponsiveVP
VPalphaMix = fracMixedSelectiveVP/fracMixedSelectiveDP

nMixedSelectiveVA = np.sum(quadrantsMixedSelective[3])
nSpeechSelectiveVA = np.sum(quadrantsMixedSelective[3]) + np.sum(quadrantsSingleSelective[3])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3]& quadrantsUseLayers[3])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3]& quadrantsUseLayers[3])
nTotalVA = quadrantTotals[3] - np.sum(~quadrantsUseLayers[3])
#fracMixedSelectiveVA = nMixedSelectiveVA/nSpeechSelectiveVA
fracMixedSelectiveVA = nMixedSelectiveVA/nSpeechResponsiveVA
VAalphaMix = fracMixedSelectiveVA/fracMixedSelectiveDP



pltMixDP = plt.Rectangle((0,0.5), 0.5, 0.5, facecolor = colorMixedSelective, alpha = DPalphaMix)
pltMixDA = plt.Rectangle((0.5,0.5), 0.5, 0.5, facecolor = colorMixedSelective, alpha = DAalphaMix)
pltMixVP = plt.Rectangle((0,0), 0.5, 0.5, facecolor = colorMixedSelective, alpha = VPalphaMix)
pltMixVA = plt.Rectangle((0.5,0), 0.5, 0.5, facecolor = colorMixedSelective, alpha = VAalphaMix)
axQuadPcts.add_artist(pltMixDP)
axQuadPcts.add_artist(pltMixDA)
axQuadPcts.add_artist(pltMixVP)
axQuadPcts.add_artist(pltMixVA)
'''
plt.annotate(f'DP\n{np.round(fracMixedSelectiveDP*100,1)}%\n n = {nSpeechSelectiveDP}', (0.25, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'w')
plt.annotate(f'DA\n{np.round(fracMixedSelectiveDA*100,1)}%\n n = {nSpeechSelectiveDA}', (0.75, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'k')
plt.annotate(f'VP\n{np.round(fracMixedSelectiveVP*100,1)}%\n n = {nSpeechSelectiveVP}', (0.25, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k')
plt.annotate(f'VA\n{np.round(fracMixedSelectiveVA*100,1)}%\n n = {nSpeechSelectiveVA}', (0.75, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k')
'''
plt.annotate(f'DP\n{np.round(fracMixedSelectiveDP*100,1)}%\nn = {nSpeechResponsiveDP}', (0.25, 0.6), fontsize=fontSizeLabels, ha='center', color='w')
plt.annotate(f'DA\n{np.round(fracMixedSelectiveDA*100,1)}%\nn = {nSpeechResponsiveDA}', (0.75, 0.6), fontsize=fontSizeLabels, ha='center', color='k')
plt.annotate(f'VP\n{np.round(fracMixedSelectiveVP*100,1)}%\nn = {nSpeechResponsiveVP}', (0.25, 0.1), fontsize=fontSizeLabels, ha='center', color='k')
plt.annotate(f'VA\n{np.round(fracMixedSelectiveVA*100,1)}%\nn = {nSpeechResponsiveVA}', (0.75, 0.1), fontsize=fontSizeLabels, ha='center', color='k')
axQuadPcts.spines["right"].set_visible(False)
axQuadPcts.spines["top"].set_visible(False)
axQuadPcts.set_aspect('equal')
plt.xticks([0, 0.5, 1], labels = np.round(quadrantBoundsAP, 2))
plt.yticks([0, 0.5, 1], labels = np.round(quadrantBoundsDV[::-1], 2))
plt.xlim([-0.1, 1.1])
plt.ylim([-0.1, 1.1])
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)


plt.sca(axByAnimal)
#fracMixedSelectiveByAnimal = quadrantsMixedSelectiveByAnimal/(quadrantsMixedSelectiveByAnimal + quadrantsSingleSelectiveByAnimal)
fracMixedSelectiveByAnimal = quadrantsMixedSelectiveByAnimal/(quadrantsSpeechResponsiveByAnimal)
fracMixedSelectiveByAnimal[quadrantsSpeechResponsiveByAnimal<=2] = np.nan
medianMixedSelectiveByAnimal = np.nanmedian(fracMixedSelectiveByAnimal, axis = 0)
quadXcoords = np.array([1,2,3,4])
plt.bar(quadXcoords, [medianMixedSelectiveByAnimal[0], medianMixedSelectiveByAnimal[1], medianMixedSelectiveByAnimal[2], medianMixedSelectiveByAnimal[3]], facecolor = colorMixedSelective)
plt.xticks(quadXcoords, ['DP', 'DA', 'VP', 'VA'])
plt.ylabel('Fraction mixed-selective')
plt.xlabel('AC regions')
for indAnimal, thisAnimal in enumerate(subjects):
    plt.plot(quadXcoords[quadrantsSpeechResponsiveByAnimal[indAnimal]>2],
             fracMixedSelectiveByAnimal[indAnimal,quadrantsSpeechResponsiveByAnimal[indAnimal]>2],
             c=colorNotSelective, marker='o', alpha=0.5)
axByAnimal.spines["right"].set_visible(False)
axByAnimal.spines["top"].set_visible(False)

titleY = 0.9

plt.sca(axMixAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0]) & useLayers)
nMixselectiveAudP = np.sum(mixedSelective[recordingAreaName == audCtxAreas[0]])
nSingleselectiveAudP = np.sum(singleSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nMixselectiveAudP, nSingleselectiveAudP, nSpeechResponsiveAudP - (nMixselectiveAudP + nSingleselectiveAudP)], colors = [colorMixedSelective, colorSingleSelective, colorNotSelective])
#plt.pie([nMixselectiveAudP, nSingleselectiveAudP], colors = [colorMixedSelective, colorSingleSelective])

axMixAudP.add_artist(circle1)
plt.title(f'AudP (n={int(nSpeechResponsiveAudP)})', y=titleY, fontsize=fontSizeLabels)
#plt.title(f'AudP,\n n = {int(nMixselectiveAudP+nSingleselectiveAudP)}', pad = 0 )


plt.sca(axMixAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1]) & useLayers)
nMixselectiveAudD = np.sum(mixedSelective[recordingAreaName == audCtxAreas[1]])
nSingleselectiveAudD = np.sum(singleSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nMixselectiveAudD, nSingleselectiveAudD, nSpeechResponsiveAudD - (nMixselectiveAudD + nSingleselectiveAudD)], colors = [colorMixedSelective, colorSingleSelective, colorNotSelective])
#plt.pie([nMixselectiveAudD, nSingleselectiveAudD], colors = [colorMixedSelective, colorSingleSelective])
axMixAudD.add_artist(circle2)
plt.title(f'AudD (n={int(nSpeechResponsiveAudD)})', y=titleY, fontsize=fontSizeLabels)
#plt.title(f'AudD,\n n = {int(nMixselectiveAudD+nSingleselectiveAudD)}', pad = 0)


plt.sca(axMixAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2]) & useLayers)
nMixselectiveAudV = np.sum(mixedSelective[recordingAreaName == audCtxAreas[2]])
nSingleselectiveAudV = np.sum(singleSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nMixselectiveAudV, nSingleselectiveAudV, nSpeechResponsiveAudV - (nMixselectiveAudV + nSingleselectiveAudV)], colors = [colorMixedSelective, colorSingleSelective, colorNotSelective])
#plt.pie([nMixselectiveAudV, nSingleselectiveAudV], colors = [colorMixedSelective, colorSingleSelective])
axMixAudV.add_artist(circle3)
plt.title(f'AudV (n={int(nSpeechResponsiveAudV)})', y=titleY, fontsize=fontSizeLabels)
#plt.title(f'AudV,\n n = {int(nMixselectiveAudV+nSingleselectiveAudV)}', pad = 0)

plt.sca(axMixTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3]) & useLayers)
nMixselectiveTeA = np.sum(mixedSelective[recordingAreaName == audCtxAreas[0]])
nSingleselectiveTeA = np.sum(singleSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nMixselectiveTeA, nSingleselectiveTeA, nSpeechResponsiveTeA - (nMixselectiveTeA + nSingleselectiveTeA)], colors = [colorMixedSelective, colorSingleSelective, colorNotSelective])
#plt.pie([nMixselectiveTeA, nSingleselectiveTeA], colors = [colorMixedSelective, colorSingleSelective])
axMixTeA.add_artist(circle4)
plt.title(f'TeA (n={int(nSpeechResponsiveTeA)})', y=titleY, fontsize=fontSizeLabels)
#plt.title(f'TeA,\n n = {int(nMixselectiveTeA+nSingleselectiveTeA)}', pad = 0)
#plt.legend(labels = ['Non-selective', 'Mixed-selective', 'Single-selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))


labelPosX = [0.01, 0.45, 0.67] # Horiz position for panel labels
labelPosY = [0.96, 0.46]    # Vert position for panel labels

axMixeSelMap.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axMixAudD.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axQuadPcts.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axByAnimal.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

plt.show()


if STATSUMMARY:
    # --ATLAS AREAS
    atlasAreasMixedSelective = np.array([nMixselectiveAudP, nMixselectiveAudD, nMixselectiveAudV, nMixselectiveTeA])
    atlasAreasSingleSeletive = np.array([nSingleselectiveAudP, nSingleselectiveAudD, nSingleselectiveAudV, nSingleselectiveTeA])
    atlasAreasSpeechResponsive = np.array([nSpeechResponsiveAudP, nSpeechResponsiveAudD, nSpeechResponsiveAudV, nSpeechResponsiveTeA])

    nCompar = 6
    ##--Test frac mixed selective
    atlasAreaComparFracMixedSelective = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(atlasAreasMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(atlasAreasSingleSeletive[indBin]), np.sum(atlasAreasSingleSeletive[x + indBin + 1])],[np.sum(atlasAreasMixedSelective[indBin]), np.sum(atlasAreasMixedSelective[x + indBin + 1])]]))
            atlasAreaComparFracMixedSelective[a] = pvalFracMixedSelective

    atlasAreaComparFracMixedSelective_speechResponsive = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(atlasAreasMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(atlasAreasSpeechResponsive[indBin])-np.sum(atlasAreasMixedSelective[indBin]), np.sum(atlasAreasSpeechResponsive[x + indBin + 1]) - np.sum(atlasAreasMixedSelective[x + indBin + 1])],[np.sum(atlasAreasMixedSelective[indBin]), np.sum(atlasAreasMixedSelective[x + indBin + 1])]]))
            atlasAreaComparFracMixedSelective_speechResponsive[a] = pvalFracMixedSelective

    atlasComparLabels = np.array(['AudP vs AudD', 'AudP vs. AudV', 'AudP vs. TeA', 'AudD vs. AudV', 'AudD vs. TeA', 'AudV vs. TeA'])

    # --QUADRANTS
    nQuadCompar = 6
    quadAlpha = 0.05/6
    quadComparLabels = np.array(['DP vs DA', 'DP vs. VP', 'DP vs. VA', 'DA vs. VP', 'DA vs. VA', 'VP vs. VA'])
    ##--Test frac mixed selective
    quadrantComparFracMixedSelective = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quadrantsSingleSelective[indBin]), np.sum(quadrantsSingleSelective[x + indBin + 1])],[np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsMixedSelective[x + indBin + 1])]]))
            quadrantComparFracMixedSelective[a] = pvalFracMixedSelective

    quadrantComparFracMixedSelective_speechResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quadrantsSpeechResponsive[indBin])-np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1]) - np.sum(quadrantsMixedSelective[x + indBin + 1])],[np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsMixedSelective[x + indBin + 1])]]))
            quadrantComparFracMixedSelective_speechResponsive[a] = pvalFracMixedSelective

    quadrantComparFracMixedSelective_soundResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quadrantsSoundResponsive[indBin])-np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsMixedSelective[x + indBin + 1])],[np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsMixedSelective[x + indBin + 1])]]))
            quadrantComparFracMixedSelective_soundResponsive[a] = pvalFracMixedSelective

    # -- Print LaTeX tables --
    reload(studyutils)
    brainAreas = ['AudP', 'AudD', 'AudV', 'TeA']
    nSpeechResponsiveEachArea = [eval('nSpeechResponsive'+brainArea) for brainArea in brainAreas]
    nMixselectiveEachArea = [eval('nMixselective'+brainArea) for brainArea in brainAreas]
    nSingleselectiveEachArea = [eval('nSingleselective'+brainArea) for brainArea in brainAreas]
    nSpeechSelectiveEachArea = np.array(nMixselectiveEachArea) + np.array(nSingleselectiveEachArea)

    pvalsAtlasAreaComparFracMixedSelective = np.full((len(brainAreas),len(brainAreas)), np.nan)
    a = -1
    for inda1, area1 in enumerate(brainAreas):
        for inda2, area2 in enumerate(brainAreas):
            if inda2 > inda1:
                a = a+1
                pvalsAtlasAreaComparFracMixedSelective[inda1,inda2] = atlasAreaComparFracMixedSelective_speechResponsive[a]

    brainRegions = ['DP', 'DA', 'VP', 'VA']
    nSpeechResponsiveEachRegion = [eval('nSpeechResponsive'+brainRegion) for brainRegion in brainRegions]
    nSoundResponsiveEachRegion = [eval('nSoundResponsive'+brainRegion) for brainRegion in brainRegions]
    nSpeechSelectiveEachRegion = [eval('nSpeechSelective'+brainRegion) for brainRegion in brainRegions]
    nMixselectiveEachRegion = [eval('nMixedSelective'+brainRegion) for brainRegion in brainRegions]


    pvalQuadrantComparFracMixedSelective = np.full((len(brainRegions),len(brainRegions)), np.nan)
    a = -1
    for inda1, area1 in enumerate(brainRegions):
        for inda2, area2 in enumerate(brainRegions):
            if inda2 > inda1:
                a = a+1
                pvalQuadrantComparFracMixedSelective[inda1,inda2] = quadrantComparFracMixedSelective_speechResponsive[a]

    table9a = studyutils.latex_table_mixselective(brainAreas, nMixselectiveEachArea, nSpeechResponsiveEachArea, nSpeechSelectiveEachArea)
    print(table9a)
    table9b = studyutils.latex_table_pvals(brainAreas, pvalsAtlasAreaComparFracMixedSelective)
    print(table9b)
    table10a = studyutils.latex_table_mixselective(brainRegions, nMixselectiveEachRegion, nSpeechResponsiveEachRegion, nSpeechSelectiveEachRegion)
    print(table10a)
    table10b = studyutils.latex_table_pvals(brainRegions, pvalQuadrantComparFracMixedSelective)
    print(table10b)



    print('--Stat Summary --')
    print('--atlas-defined areas--')
    print(f'AudP: n mixed selective = {nMixselectiveAudP} ({np.round(nMixselectiveAudP/(nMixselectiveAudP+nSingleselectiveAudP)*100, 1)}% of speech selective, {np.round(nMixselectiveAudP/(nSpeechResponsiveAudP)*100, 1)}% of speech responsive)')
    print(f'AudD: n mixed selective = {nMixselectiveAudD} ({np.round(nMixselectiveAudD/(nMixselectiveAudD+nSingleselectiveAudD)*100, 1)}% of speech selective, {np.round(nMixselectiveAudD/(nSpeechResponsiveAudD)*100, 1)}% of speech responsive)')
    print(f'AudV: n mixed selective = {nMixselectiveAudV} ({np.round(nMixselectiveAudV/(nMixselectiveAudV+nSingleselectiveAudV)*100, 1)}% of speech selective, {np.round(nMixselectiveAudV/(nSpeechResponsiveAudV)*100, 1)}% of speech responsive)')
    print(f'TeA: n mixed selective = {nMixselectiveTeA} ({np.round(nMixselectiveTeA/(nMixselectiveTeA+nSingleselectiveTeA)*100, 1)}% of speech selective, {np.round(nMixselectiveTeA/(nSpeechResponsiveTeA)*100, 1)}% of speech responsive)')
    print('--Frac mixed selective of speech selective cells --')
    print(f'AudP vs AudD p = {np.round(atlasAreaComparFracMixedSelective[0],3)}')
    print(f'AudP vs AudV p = {np.round(atlasAreaComparFracMixedSelective[1],3)}')
    print(f'AudP vs TeA p = {np.round(atlasAreaComparFracMixedSelective[2],3)}')
    print(f'AudD vs AudV p = {np.round(atlasAreaComparFracMixedSelective[3],3)}')
    print(f'AudD vs TeA p = {np.round(atlasAreaComparFracMixedSelective[4],3)}')
    print(f'AudV vs TeA p = {np.round(atlasAreaComparFracMixedSelective[5],3)}')
    print('--Frac mixed selective of speech responsive --')
    print(f'AudP vs AudD p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[0],3)}')
    print(f'AudP vs AudV p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[1],3)}')
    print(f'AudP vs TeA p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[2],3)}')
    print(f'AudD vs AudV p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[3],3)}')
    print(f'AudD vs TeA p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[4],3)}')
    print(f'AudV vs TeA p = {np.round(atlasAreaComparFracMixedSelective_speechResponsive[5],3)}')

    print('--quadrant-defined areas--')
    print(f'DP Quadrant: n mixed selective = {np.sum(quadrantsMixedSelective[0])} ({np.round((np.sum(quadrantsMixedSelective[0])/(np.sum(quadrantsMixedSelective[0]) + np.sum(quadrantsSingleSelective[0])))*100,1)}% of speech selective, {np.round((np.sum(quadrantsMixedSelective[0])/np.sum(quadrantsSpeechResponsive[0]))*100,1)}% of speech responsive)')
    print(f'DA Quadrant: n mixed selective = {np.sum(quadrantsMixedSelective[1])} ({np.round((np.sum(quadrantsMixedSelective[1])/(np.sum(quadrantsMixedSelective[1]) + np.sum(quadrantsSingleSelective[1])))*100,1)}% of speech selective, {np.round((np.sum(quadrantsMixedSelective[1])/np.sum(quadrantsSpeechResponsive[1]))*100,1)}% of speech responsive)')
    print(f'VP Quadrant:  n mixed selective = {np.sum(quadrantsMixedSelective[2])} ({np.round((np.sum(quadrantsMixedSelective[2])/(np.sum(quadrantsMixedSelective[2]) + np.sum(quadrantsSingleSelective[2])))*100,1)}% of speech selective, {np.round((np.sum(quadrantsMixedSelective[3])/np.sum(quadrantsSpeechResponsive[2]))*100,1)}% of speech responsive)')
    print(f'VA Quadrant: n mixed selective = {np.sum(quadrantsMixedSelective[3])} ({np.round((np.sum(quadrantsMixedSelective[2])/(np.sum(quadrantsMixedSelective[3]) + np.sum(quadrantsSingleSelective[3])))*100,1)}% of speech selective, {np.round((np.sum(quadrantsMixedSelective[3])/np.sum(quadrantsSpeechResponsive[3]))*100,1)}% of speech responsive)')
    print('--Frac mixed selective of speech selective cells --')
    print(f'DP vs DA p = {np.round(quadrantComparFracMixedSelective[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracMixedSelective[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracMixedSelective[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracMixedSelective[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracMixedSelective[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracMixedSelective[5],3)}')
    print('--Frac mixed selective of speech responsive --')
    print(f'DP vs DA p = {np.round(quadrantComparFracMixedSelective_speechResponsive[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracMixedSelective_speechResponsive[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracMixedSelective_speechResponsive[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracMixedSelective_speechResponsive[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracMixedSelective_speechResponsive[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracMixedSelective_speechResponsive[5],3)}')

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
