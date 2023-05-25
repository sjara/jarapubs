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

reload(figparams)

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
shuffledDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_mixed_selectivity_wResponsive' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 5.25] # In inches
STATSUMMARY = 0

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
fontSizePanel = figparams.fontSizePanel
colorMap = cm.get_cmap('Greens')
newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
colorMixedSelective = cp.TangoPalette['ScarletRed2']
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
excludeCells = figData['excludeCells']
bestSelectivityIndexVot = figData['bestSelectivityIndexVot']
bestSelectivityIndexFt = figData['bestSelectivityIndexFt']
pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
shuffledVotBest = figData['shuffledVotBest']
shuffledFtBest = figData['shuffledFtBest']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
avgShuffledSIVot = np.mean(shuffledVotBest[~excludeCells & speechResponsive],1)
avgShuffledSIFt = np.mean(shuffledFtBest[~excludeCells & speechResponsive], 1)


plt.figure()
gsMain = gridspec.GridSpec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axMixeSelMap = plt.subplot(gsMain[:,0])
axMixedBins = gsMain[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45])
axMixDV = plt.subplot(axMixedBins[0,0])
axMixAP = plt.subplot(axMixedBins[1,0])
axMixDonuts = gsMain[0:,1].subgridspec(4, 1, hspace = 0.4)
axMixAudP = plt.subplot(axMixDonuts[1,0])
axMixAudD = plt.subplot(axMixDonuts[0,0])
axMixAudV = plt.subplot(axMixDonuts[2,0])
axMixTeA = plt.subplot(axMixDonuts[3,0])

gsMain.update(left=0.08, right=0.96, top=0.92, bottom=0.1, wspace=0.25, hspace=0.4)
plt.subplots_adjust(top = 0.9, bottom = 0.1, hspace = 0.45, left = 0.05)

nBins = 8
if nBins == 10:
    nCompar = 45
elif nBins == 9:
    nCompar = 36
elif nBins == 8:
    nCompar = 28

binSizeDV = (np.max(y_coords[speechResponsive & ~excludeCells]) - np.min(y_coords[speechResponsive & ~excludeCells]))/nBins
binsDV = np.arange(np.min(y_coords[speechResponsive & ~excludeCells]), np.max(y_coords[speechResponsive & ~excludeCells]), binSizeDV)

binSizeAP = (np.max(z_coords[speechResponsive & ~excludeCells]) - np.min(z_coords[speechResponsive & ~excludeCells]))/nBins
binsAP = np.arange(np.min(z_coords[speechResponsive & ~excludeCells]), np.max(z_coords[speechResponsive & ~excludeCells]), binSizeAP)

votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive & ~excludeCells
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive & ~excludeCells
singleSelective = np.logical_xor(votSelective, ftSelective)
mixedSelective = votSelective & ftSelective

quantilesVOT_DV = []
quantilesFT_DV = []
quantilesVOT_AP = []
quantilesFT_AP = []
quantilesVotSelective_DV = []
quantilesVotSelective_AP = []
quantilesFtSelective_DV = []
quantilesFtSelective_AP = []
quantilesSingleSelective_AP = []
quantilesSingleSelective_DV = []
quantilesMixedSelective_AP = []
quantilesMixedSelective_DV = []
quantilesResponsive_AP = []
quantilesResponsive_DV = []

for indBin, thisBin in enumerate(binsDV):
    if indBin < len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords < binsDV[indBin+1])
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords < binsAP[indBin+1])
    elif indBin == len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords <= np.max(y_coords[speechResponsive & ~excludeCells]))
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords <= np.max(z_coords[speechResponsive & ~excludeCells]))
    quantilesVOT_DV.append(bestSelectivityIndexVot[thisQuantileDV & speechResponsive & ~excludeCells])
    quantilesFT_DV.append(bestSelectivityIndexFt[thisQuantileDV & speechResponsive & ~excludeCells])
    quantilesVOT_AP.append(bestSelectivityIndexVot[thisQuantileAP & speechResponsive & ~excludeCells])
    quantilesFT_AP.append(bestSelectivityIndexFt[thisQuantileAP & speechResponsive & ~excludeCells])
    quantilesVotSelective_AP.append(votSelective[thisQuantileAP])
    quantilesFtSelective_AP.append(ftSelective[thisQuantileAP])
    quantilesVotSelective_DV.append(votSelective[thisQuantileDV])
    quantilesFtSelective_DV.append(ftSelective[thisQuantileDV])
    quantilesSingleSelective_AP.append(singleSelective[thisQuantileAP])
    quantilesSingleSelective_DV.append(singleSelective[thisQuantileDV])
    quantilesMixedSelective_AP.append(mixedSelective[thisQuantileAP])
    quantilesMixedSelective_DV.append(mixedSelective[thisQuantileDV])
    quantilesResponsive_AP.append(speechResponsive[thisQuantileAP])
    quantilesResponsive_DV.append(speechResponsive[thisQuantileDV])


APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)

plt.sca(axMixeSelMap)
nonSel = plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells & ~mixedSelective & ~singleSelective], y_coords[speechResponsive & ~excludeCells & ~mixedSelective & ~singleSelective], c = colorNotSelective, s = 6)
singSel = plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells & singleSelective], y_coords[speechResponsive & ~excludeCells & singleSelective], c = colorSingleSelective, s = 6)
mixSel = plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells & mixedSelective], y_coords[speechResponsive & ~excludeCells & mixedSelective], c = colorMixedSelective, s = 6)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.legend(['Non-selective', 'Single-selective', 'Mixed-selective'], loc = "upper right", markerscale = 3, bbox_to_anchor = (1.2, 1.1))
#plt.title('Mixed selectivity', fontsize = fontSizeTitles)
axMixeSelMap.spines["right"].set_visible(False)
axMixeSelMap.spines["top"].set_visible(False)

binCountsMixed_AP = np.zeros(nBins)
binCountsMixed_DV = np.zeros(nBins)
binCountsSingle_AP = np.zeros(nBins)
binCountsSingle_DV = np.zeros(nBins)
binCountsNonSel_AP = np.zeros(nBins)
binCountsNonSel_DV = np.zeros(nBins)


# Not showing non-responsive
for indBin, thisBin in enumerate(quantilesMixedSelective_AP):
    nSelectiveAP = np.sum(quantilesMixedSelective_AP[indBin]) + np.sum(quantilesSingleSelective_AP[indBin])
    nSelectiveDV = np.sum(quantilesSingleSelective_DV[indBin]) + np.sum(quantilesMixedSelective_DV[indBin])
    binCountsMixed_AP[indBin] = np.sum(quantilesMixedSelective_AP[indBin])/nSelectiveAP
    binCountsMixed_DV[indBin] = np.sum(quantilesMixedSelective_DV[indBin])/nSelectiveDV
    binCountsSingle_AP[indBin] = np.sum(quantilesSingleSelective_AP[indBin])/nSelectiveAP
    binCountsSingle_DV[indBin] = np.sum(quantilesSingleSelective_DV[indBin])/nSelectiveDV

binCounts_AP = {"Mixed": binCountsMixed_AP, "Single": binCountsSingle_AP}
binCounts_DV = {"Mixed": binCountsMixed_DV, "Single": binCountsSingle_DV}

'''

# Showing non-responsive
for indBin, thisBin in enumerate(quantilesMixedSelective_AP):
    nSelectiveAP = np.sum(quantilesMixedSelective_AP[indBin]) + np.sum(quantilesSingleSelective_AP[indBin])
    nSelectiveDV = np.sum(quantilesSingleSelective_DV[indBin]) + np.sum(quantilesMixedSelective_DV[indBin])
    nResponsiveAP = np.sum(quantilesResponsive_AP[indBin])
    nResponsiveDV = np.sum(quantilesResponsive_DV[indBin])
    binCountsMixed_AP[indBin] = np.sum(quantilesMixedSelective_AP[indBin])/nResponsiveAP
    binCountsMixed_DV[indBin] = np.sum(quantilesMixedSelective_DV[indBin])/nResponsiveDV
    binCountsSingle_AP[indBin] = np.sum(quantilesSingleSelective_AP[indBin])/nResponsiveAP
    binCountsSingle_DV[indBin] = np.sum(quantilesSingleSelective_DV[indBin])/nResponsiveDV
    binCountsNonSel_AP[indBin] =  (nResponsiveAP-nSelectiveAP)/nResponsiveAP
    binCountsNonSel_DV[indBin] =  (nResponsiveDV-nSelectiveDV)/nResponsiveDV

binCounts_AP = {"Mixed": binCountsMixed_AP, "Single": binCountsSingle_AP, "Non-selective": binCountsNonSel_AP}
binCounts_DV = {"Mixed": binCountsMixed_DV, "Single": binCountsSingle_DV, "Non-selective": binCountsNonSel_DV}
'''

plt.sca(axMixDV)
mixedBarsDV = plt.barh(binsDV, binCountsMixed_DV, height = binSizeDV/1.75, color = colorMixedSelective) #, color = [colorMixedSelective, colorSingleSelective]
singleBarsDV = plt.barh(binsDV, binCountsSingle_DV, height = binSizeDV/1.75, left = binCountsMixed_DV, color = colorSingleSelective)
#nonSelBarsDV = plt.barh(binsDV, binCountsNonSel_DV, height = binSizeDV/1.75, left = (binCountsMixed_DV+binCountsSingle_DV), color = colorNotSelective)
plt.ylim(190, 40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Fraction of Selective Cells', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
axMixDV.spines["right"].set_visible(False)
axMixDV.spines["top"].set_visible(False)


plt.sca(axMixAP)
mixedBarsAP = plt.bar(binsAP, binCountsMixed_AP, width = binSizeAP/1.75, color = colorMixedSelective)
singleBarsAP = plt.bar(binsAP, binCountsSingle_AP, width = binSizeAP/1.75, bottom = binCountsMixed_AP, color = colorSingleSelective)
#nonSelBarsAP = plt.bar(binsAP, binCountsNonSel_AP, width = binSizeAP/1.75, bottom = (binCountsMixed_AP + binCountsSingle_AP), color = colorNotSelective)
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.ylabel('Fraction of Selective Cells', fontsize = fontSizeLabels)
axMixAP.spines["right"].set_visible(False)
axMixAP.spines["top"].set_visible(False)


plt.sca(axMixAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & ~excludeCells & (recordingAreaName == audCtxAreas[0] ))
nMixselectiveAudP = np.sum(mixedSelective[recordingAreaName == audCtxAreas[0]])
nSingleselectiveAudP = np.sum(singleSelective[recordingAreaName == audCtxAreas[0]])
#plt.pie([nSpeechResponsiveAudP - (nMixselectiveAudP + nSingleselectiveAudP),  nMixselectiveAudP, nSingleselectiveAudP], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudP, nSingleselectiveAudP], colors = [colorMixedSelective, colorSingleSelective])

axMixAudP.add_artist(circle1)
#plt.title(f'AudP,\n n = {int(nSpeechResponsiveAudP)}')
plt.title(f'AudP,\n n = {int(nMixselectiveAudP+nSingleselectiveAudP)}', pad = 0 )


plt.sca(axMixAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & ~excludeCells & (recordingAreaName == audCtxAreas[1] ))
nMixselectiveAudD = np.sum(mixedSelective[recordingAreaName == audCtxAreas[1]])
nSingleselectiveAudD = np.sum(singleSelective[recordingAreaName == audCtxAreas[1]])
#plt.pie([nSpeechResponsiveAudD - (nMixselectiveAudD + nSingleselectiveAudD),  nMixselectiveAudD, nSingleselectiveAudD], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudD, nSingleselectiveAudD], colors = [colorMixedSelective, colorSingleSelective])
axMixAudD.add_artist(circle2)
#plt.title(f'AudD,\n n = {int(nSpeechResponsiveAudD)}')
plt.title(f'AudD,\n n = {int(nMixselectiveAudD+nSingleselectiveAudD)}', pad = 0)


plt.sca(axMixAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & ~excludeCells & (recordingAreaName == audCtxAreas[2] ))
nMixselectiveAudV = np.sum(mixedSelective[recordingAreaName == audCtxAreas[2]])
nSingleselectiveAudV = np.sum(singleSelective[recordingAreaName == audCtxAreas[2]])
#plt.pie([nSpeechResponsiveAudV - (nMixselectiveAudV + nSingleselectiveAudV),  nMixselectiveAudV, nSingleselectiveAudV], colors = [colorNotSelective, colorMixedSelective, colorSingleSelective])
plt.pie([nMixselectiveAudV, nSingleselectiveAudV], colors = [colorMixedSelective, colorSingleSelective])
axMixAudV.add_artist(circle3)
#plt.title(f'AudV,\n n = {int(nSpeechResponsiveAudV)}')
plt.title(f'AudV,\n n = {int(nMixselectiveAudV+nSingleselectiveAudV)}', pad = 0)

plt.sca(axMixTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & ~excludeCells & (recordingAreaName == audCtxAreas[3] ))
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
axMixDV.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axMixAP.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
