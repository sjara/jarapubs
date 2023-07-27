"""
This creates figure 1 for 2022paspeech:
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
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import celldatabase
from jaratoolbox import colorpalette as cp
import matplotlib.colors as colors
import matplotlib.cm as cm
from scipy import stats
import studyparams
import figparams
from importlib import reload
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import findfont, FontProperties
font = findfont(FontProperties(family = ['Helvetica']))

reload(figparams)

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
shuffledDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_mixed_selectivity' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [7.5, 5.25] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
fontSizePanel = figparams.fontSizePanel
colorMap = cm.get_cmap('Greens')
newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
colorMixedSelective = cp.TangoPalette['Orange2']
colorSingleSelective = cp.TangoPalette['SkyBlue2']
colorNotSelective = cp.TangoPalette['Aluminium3']


figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
z_coords_jittered = figData['z_coords_jittered']
x_coords_jittered = figData['x_coords_jittered']
speechResponsive = figData['speechResponsive']
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
quadrantBoundsDV_AAtransform = figData['quadrantBoundsDV_AAtransform']
quadrantBoundsAP_AAtransform = figData['quadrantBoundsAP_AAtransform']
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


plt.figure()
gsMain = gridspec.GridSpec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axMixeSelMap = plt.subplot(gsMain[:,0])
axQuadSummary = gsMain[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45])
axQuadPcts = plt.subplot(axQuadSummary[0,0])
axByAnimal = plt.subplot(axQuadSummary[1,0])
axMixDonuts = gsMain[0:,1].subgridspec(4, 1, hspace = 0.4)
axMixAudP = plt.subplot(axMixDonuts[1,0])
axMixAudD = plt.subplot(axMixDonuts[0,0])
axMixAudV = plt.subplot(axMixDonuts[2,0])
axMixTeA = plt.subplot(axMixDonuts[3,0])

gsMain.update(left=0.08, right=0.96, top=0.92, bottom=0.1, wspace=0.25, hspace=0.4)
plt.subplots_adjust(top = 0.9, bottom = 0.1, hspace = 0.45, left = 0.05)


votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive
singleSelective = np.logical_xor(votSelective, ftSelective)
mixedSelective = votSelective & ftSelective



APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)

plt.sca(axMixeSelMap)
#nonSel = plt.scatter(z_coords_jittered[speechResponsive & ~mixedSelective & ~singleSelective], y_coords[speechResponsive & ~mixedSelective & ~singleSelective], c = colorNotSelective, s = 6)
singSel = plt.scatter(z_coords_jittered[speechResponsive & singleSelective], y_coords[speechResponsive & singleSelective], c = colorSingleSelective, s = 6)
mixSel = plt.scatter(z_coords_jittered[speechResponsive & mixedSelective], y_coords[speechResponsive & mixedSelective], c = colorMixedSelective, s = 6)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.legend(handles = [singSel, mixSel], labels = ['Single-selective', 'Mixed-selective'], loc = "upper right", markerscale = 3, bbox_to_anchor = (1.2, 1.1))
#plt.title('Mixed selectivity', fontsize = fontSizeTitles)
axMixeSelMap.spines["right"].set_visible(False)
axMixeSelMap.spines["top"].set_visible(False)



plt.sca(axQuadPcts)
nMixedSelectiveDP = np.sum(quadrantsMixedSelective[0])
nSpeechSelectiveDP = np.sum(quadrantsMixedSelective[0]) + np.sum(quadrantsSingleSelective[0])
nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0])
nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0])
nTotalDP = quadrantTotals[0]
fracMixedSelectiveDP = nMixedSelectiveDP/nSpeechSelectiveDP
DPalphaMix = 1

nMixedSelectiveDA = np.sum(quadrantsMixedSelective[1])
nSpeechSelectiveDA = np.sum(quadrantsMixedSelective[1]) + np.sum(quadrantsSingleSelective[1])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1])
nTotalDA = quadrantTotals[1]
fracMixedSelectiveDA = nMixedSelectiveDA/nSpeechSelectiveDA
DAalphaMix = fracMixedSelectiveDA/fracMixedSelectiveDP

nMixedSelectiveVP = np.sum(quadrantsMixedSelective[2])
nSpeechSelectiveVP = np.sum(quadrantsMixedSelective[2]) + np.sum(quadrantsSingleSelective[2])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2])
nTotalVP = quadrantTotals[2]
fracMixedSelectiveVP = nMixedSelectiveVP/nSpeechSelectiveVP
VPalphaMix = fracMixedSelectiveVP/fracMixedSelectiveDP

nMixedSelectiveVA = np.sum(quadrantsMixedSelective[3])
nSpeechSelectiveVA = np.sum(quadrantsMixedSelective[3]) + np.sum(quadrantsSingleSelective[3])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3])
nTotalVA = quadrantTotals[3]
fracMixedSelectiveVA = nMixedSelectiveVA/nSpeechSelectiveVA
VAalphaMix = fracMixedSelectiveVA/fracMixedSelectiveDP



pltMixDP = plt.Rectangle((0,0.5), 0.5, 0.5, facecolor = colorMixedSelective, alpha = DPalphaMix)
pltMixDA = plt.Rectangle((0.5,0.5), 0.5, 0.5, facecolor = colorMixedSelective, alpha = DAalphaMix)
pltMixVP = plt.Rectangle((0,0), 0.5, 0.5, facecolor = colorMixedSelective, alpha = VPalphaMix)
pltMixVA = plt.Rectangle((0.5,0), 0.5, 0.5, facecolor = colorMixedSelective, alpha = VAalphaMix)
axQuadPcts.add_artist(pltMixDP)
axQuadPcts.add_artist(pltMixDA)
axQuadPcts.add_artist(pltMixVP)
axQuadPcts.add_artist(pltMixVA)
plt.annotate(f'{np.round(fracMixedSelectiveDP*100,1)}% DP\n n = {nSpeechSelectiveDP}', (0.1, 0.65), fontsize = fontSizeTicks, color = 'w')
plt.annotate(f'{np.round(fracMixedSelectiveDA*100,1)}% DA\n n = {nSpeechSelectiveDA}', (0.6, 0.65), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracMixedSelectiveVP*100,1)}% VP\n n = {nSpeechSelectiveVP}', (0.1, 0.15), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracMixedSelectiveVA*100,1)}% VA\n n = {nSpeechSelectiveVA}', (0.6, 0.15), fontsize = fontSizeTicks, color = 'k')
plt.axis('off')



plt.sca(axByAnimal)
plt.bar([1, 2, 3, 4], [fracMixedSelectiveDP, fracMixedSelectiveDA, fracMixedSelectiveVA, fracMixedSelectiveVP], facecolor = colorMixedSelective)
plt.xticks([1,2,3,4], ['DP', 'DA', 'VA', 'VP'])
plt.ylabel('Fraction Mixed selective')
fracMixedSelectiveByAnimal = quadrantsMixedSelectiveByAnimal/(quadrantsMixedSelectiveByAnimal + quadrantsSingleSelectiveByAnimal)
plt.plot([1,2], [fracMixedSelectiveByAnimal[:,0], fracMixedSelectiveByAnimal[:,1]], c = colorNotSelective, alpha = 0.5)
plt.plot([1,3], [fracMixedSelectiveByAnimal[:,0], fracMixedSelectiveByAnimal[:,3]], c = colorNotSelective, alpha = 0.5)
plt.plot([1,4], [fracMixedSelectiveByAnimal[:,0], fracMixedSelectiveByAnimal[:,2]], c = colorNotSelective, alpha = 0.5)
axByAnimal.spines["right"].set_visible(False)
axByAnimal.spines["top"].set_visible(False)


plt.sca(axMixAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0] ))
nMixselectiveAudP = np.sum(mixedSelective[recordingAreaName == audCtxAreas[0]])
nSingleselectiveAudP = np.sum(singleSelective[recordingAreaName == audCtxAreas[0]])
#plt.pie([nSpeechResponsiveAudP - (nMixselectiveAudP + nSingleselectiveAudP),  nMixselectiveAudP, nSingleselectiveAudP], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudP, nSingleselectiveAudP], colors = [colorMixedSelective, colorSingleSelective])

axMixAudP.add_artist(circle1)
#plt.title(f'AudP,\n n = {int(nSpeechResponsiveAudP)}')
plt.title(f'AudP,\n n = {int(nMixselectiveAudP+nSingleselectiveAudP)}', pad = 0 )


plt.sca(axMixAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1] ))
nMixselectiveAudD = np.sum(mixedSelective[recordingAreaName == audCtxAreas[1]])
nSingleselectiveAudD = np.sum(singleSelective[recordingAreaName == audCtxAreas[1]])
#plt.pie([nSpeechResponsiveAudD - (nMixselectiveAudD + nSingleselectiveAudD),  nMixselectiveAudD, nSingleselectiveAudD], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudD, nSingleselectiveAudD], colors = [colorMixedSelective, colorSingleSelective])
axMixAudD.add_artist(circle2)
#plt.title(f'AudD,\n n = {int(nSpeechResponsiveAudD)}')
plt.title(f'AudD,\n n = {int(nMixselectiveAudD+nSingleselectiveAudD)}', pad = 0)


plt.sca(axMixAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2] ))
nMixselectiveAudV = np.sum(mixedSelective[recordingAreaName == audCtxAreas[2]])
nSingleselectiveAudV = np.sum(singleSelective[recordingAreaName == audCtxAreas[2]])
#plt.pie([nSpeechResponsiveAudV - (nMixselectiveAudV + nSingleselectiveAudV),  nMixselectiveAudV, nSingleselectiveAudV], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudV, nSingleselectiveAudV], colors = [colorMixedSelective, colorSingleSelective])
axMixAudV.add_artist(circle3)
#plt.title(f'AudV,\n n = {int(nSpeechResponsiveAudV)}')
plt.title(f'AudV,\n n = {int(nMixselectiveAudV+nSingleselectiveAudV)}', pad = 0)

plt.sca(axMixTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3] ))
nMixselectiveTeA = np.sum(mixedSelective[recordingAreaName == audCtxAreas[0]])
nSingleselectiveTeA = np.sum(singleSelective[recordingAreaName == audCtxAreas[0]])
#plt.pie([nSpeechResponsiveTeA - (nMixselectiveTeA + nSingleselectiveTeA),  nMixselectiveTeA, nSingleselectiveTeA], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveTeA, nSingleselectiveTeA], colors = [colorMixedSelective, colorSingleSelective])
axMixTeA.add_artist(circle4)
plt.title(f'TeA,\n n = {int(nMixselectiveTeA+nSingleselectiveTeA)}', pad = 0)
#plt.title(f'TeA,\n n = {int(nSpeechResponsiveTeA)}')
#plt.legend(labels = ['Non-selective', 'Mixed-selective', 'Single-selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))


labelPosX = [0.05, 0.45, 0.63] # Horiz position for panel labels
labelPosY = [0.94, 0.46]    # Vert position for panel labels

axMixeSelMap.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axMixAudD.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axQuadPcts.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axByAnimal.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)