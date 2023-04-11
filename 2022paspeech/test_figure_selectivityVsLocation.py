"""
This creates a figure depicting the selectivity to speech features and how it is organized in DV - AP space. 
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
figFilename = 'figure_selectivityIndices' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [12, 7.5] # In inches


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
colorMap = cm.get_cmap('Greens')
newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
colorPts = cp.TangoPalette['SkyBlue1']
colorRandSI = cp.TangoPalette['Aluminium4']



figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath)
shuffledDataFullPath = os.path.join(figDataDir,shuffledDataFile)
shuffledData = np.load(shuffledDataFullPath, allow_pickle=True)

shuffledSI = shuffledData['shuffledSI']
shuffledSiEachCell = np.concatenate((shuffledSI[0], shuffledSI[1], shuffledSI[2], shuffledSI[3], shuffledSI[4], shuffledSI[5], shuffledSI[6]))
avgShuffledSiEachCell = np.mean(shuffledSiEachCell, 1)
avgShuffledSI = np.nanmean(avgShuffledSiEachCell)


x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
z_coords_jittered = figData['z_coords_jittered']
x_coords_jittered = figData['x_coords_jittered']
speechResponsive = figData['speechResponsive']
excludeCells = figData['excludeCells']
bestSelectivityIndexVot = figData['bestSelectivityIndexVot']
bestSelectivityIndexFt = figData['bestSelectivityIndexFt']
#meanRandSI = np.mean(figData['randSI'])

pvalPermutationtestFt = np.ones(len(bestSelectivityIndexFt))
pvalPermutationtestVot = np.ones(len(bestSelectivityIndexFt))

for indCell, thisCell in enumerate(bestSelectivityIndexFt):
    pvalPermutationtestFt[indCell] = np.mean(shuffledSiEachCell[indCell,:] >= bestSelectivityIndexFt[indCell])
    pvalPermutationtestVot[indCell] = np.mean(shuffledSiEachCell[indCell,:] >= bestSelectivityIndexVot[indCell])


plt.figure()
gsMain = gridspec.GridSpec(2, 1)

gsVOT = gsMain[0].subgridspec(1, 4, width_ratios = [0.26, 0.08, 0.33, 0.33])
axColorMapVOT = plt.subplot(gsVOT[0,0])
axVotDV = plt.subplot(gsVOT[0,2])
axVotAP = plt.subplot(gsVOT[0,3])

gsFT = gsMain[1].subgridspec(1, 4, width_ratios = [0.26, 0.08, 0.33, 0.33])
axColorMapFT = plt.subplot(gsFT[0,0])
axFtDV = plt.subplot(gsFT[0,2])
axFtAP = plt.subplot(gsFT[0,3])
gsMain.update(left=0.05, right=0.98, top=0.9, bottom=0.1, wspace=0.6, hspace=0.3)

nBins = 10
binSizeDV = (np.max(y_coords[speechResponsive & ~excludeCells]) - np.min(y_coords[speechResponsive & ~excludeCells]))/nBins
binsDV = np.arange(np.min(y_coords[speechResponsive & ~excludeCells]), np.max(y_coords[speechResponsive & ~excludeCells]), binSizeDV)

binSizeAP = (np.max(z_coords[speechResponsive & ~excludeCells]) - np.min(z_coords[speechResponsive & ~excludeCells]))/nBins
binsAP = np.arange(np.min(z_coords[speechResponsive & ~excludeCells]), np.max(z_coords[speechResponsive & ~excludeCells]), binSizeAP)

quantilesVOT_DV = []
quantilesFT_DV = []
quantilesVOT_AP = []
quantilesFT_AP = []

for indBin, thisBin in enumerate(binsDV):
    if indBin < len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords < binsDV[indBin+1])
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords < binsAP[indBin+1])
    elif indBin == len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords < np.max(y_coords[speechResponsive & ~excludeCells]))
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords < np.max(z_coords[speechResponsive & ~excludeCells]))
    quantilesVOT_DV.append(bestSelectivityIndexVot[thisQuantileDV & speechResponsive & ~excludeCells])
    quantilesFT_DV.append(bestSelectivityIndexFt[thisQuantileDV & speechResponsive & ~excludeCells])
    quantilesVOT_AP.append(bestSelectivityIndexVot[thisQuantileAP & speechResponsive & ~excludeCells])
    quantilesFT_AP.append(bestSelectivityIndexFt[thisQuantileAP & speechResponsive & ~excludeCells])


APtickLocs = np.array([176, 196, 216])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([190, 140, 90])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
#plt.suptitle('VOT selectivity by location', fontsize = fontSizeTitles)

plt.sca(axColorMapVOT)
# np.random.randn(np.sum(a1RateCoder))*0.1+1
plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells], y_coords[speechResponsive & ~excludeCells], c = bestSelectivityIndexVot[speechResponsive & ~excludeCells], cmap = newMap, s = 3)
plt.ylim(215,60)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Dorsal - Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Anterior - Posterior (mm)', fontsize = fontSizeLabels)
plt.colorbar(label = 'VOT Selectivity Index', location = 'bottom', pad = 0.2, shrink = 0.8)
plt.title('VOT selectivity map', fontsize = fontSizeTitles)
axColorMapVOT.spines["right"].set_visible(False)
axColorMapVOT.spines["top"].set_visible(False)
#plt.sca(axHistVOT_DV)
#plt.annotate('A', xy= (160, -0.75), xycoords = ('data', 'axes fraction'), xytext = (170 ,-0.75), textcoords = ('data', 'axes fraction'), arrowprops=dict(arrowstyle="<-") )
#plt.annotate('D', xy = (160, -0.75), xycoords = ('data', 'axes fraction'),  xytext = (160 ,-0.55), textcoords = ('data', 'axes fraction'), arrowprops=dict(arrowstyle="<-") )
plt.annotate('D', xy = (170, -0.75), xycoords = ('data', 'axes fraction'),  xytext = (160 ,-0.55), textcoords = ('data', 'axes fraction'), arrowprops=dict(arrowstyle="<->", connectionstyle = "angle, angleA = -90, angleB = 0, rad=0", lw = 2), fontsize = fontSizeLabels, fontweight = 'bold')
plt.annotate('A', xy =(170 ,-0.77), xycoords = ('data', 'axes fraction'), fontsize = fontSizeLabels, fontweight = 'bold')
'''
plt.annotate('V', xy = (153, -0.75), xycoords = ('data', 'axes fraction'),  xytext = (160 ,-0.95), textcoords = ('data', 'axes fraction'), arrowprops=dict(arrowstyle="<->", connectionstyle = "angle, angleA = 90, angleB = 180, rad=0", lw = 2), fontsize = fontSizeLabels, fontweight = 'bold')
plt.annotate('P', xy =(150 ,-0.77), xycoords = ('data', 'axes fraction'), fontsize = fontSizeLabels, fontweight = 'bold')
'''
#plt.text(170,100, "AudD\nAudP", size=15, bbox = dict(boxstyle="rarrow,pad=0.3", fc=cp.TangoPalette['Aluminium2'], ec=cp.TangoPalette['Aluminium3'], lw=2))
#plt.text(170,130, "AudP\nAudV", size=15, bbox = dict(boxstyle="rarrow,pad=0.3", fc=cp.TangoPalette['Aluminium2'], ec=cp.TangoPalette['Aluminium3'], lw=2))
#plt.text(170,150, "AudV\nTea", size=15, bbox = dict(boxstyle="rarrow,pad=0.3", fc=cp.TangoPalette['Aluminium2'], ec=cp.TangoPalette['Aluminium3'], lw=2))



plt.sca(axVotDV)

plt.scatter(bestSelectivityIndexVot[speechResponsive & ~excludeCells], y_coords[speechResponsive & ~excludeCells], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesVOT_DV, positions = binsDV, vert = False, widths = 6, boxprops = dict(linewidth = 2), medianprops = dict(linewidth = 2, color = 'k'), showfliers = False)
plt.plot([avgShuffledSI, avgShuffledSI] ,[215,60], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.ylim(215, 60)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('VOT Selectivity Index', fontsize = fontSizeLabels)
plt.ylabel('Dorsal - Ventral (mm)', fontsize = fontSizeLabels)
axVotDV.spines["right"].set_visible(False)
axVotDV.spines["top"].set_visible(False)




plt.sca(axVotAP)
plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells], bestSelectivityIndexVot[speechResponsive & ~excludeCells], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesVOT_AP, positions = binsAP, vert = True, widths = 2.5, boxprops = dict(linewidth = 2), medianprops = dict(linewidth = 2, color = 'k'), showfliers = False)
plt.plot([165,225], [avgShuffledSI, avgShuffledSI], linestyle = '--', linewidth = 2, c =colorRandSI)
#plt.ylim(215, 60)
#plt.yticks([])
#ax2.invert_yaxis()
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.xlabel('Anterior - Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('VOT Selectivity Index', fontsize = fontSizeLabels)
axVotAP.spines["right"].set_visible(False)
axVotAP.spines["top"].set_visible(False)



#plt.suptitle('FT selectivity by location', fontsize = fontSizeTitles)
& (pvalPermutationtestFt<0.05)
plt.sca(axColorMapFT)
plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells], y_coords[speechResponsive & ~excludeCells], c = bestSelectivityIndexFt[speechResponsive & ~excludeCells], cmap = newMap, s = 3)
plt.ylim(215,60)
plt.xlim(165,225)
plt.ylabel('Dorsal - Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Anterior - Posterior (mm)', fontsize = fontSizeLabels)
plt.xticks(APtickLocs, APtickLabels)
plt.yticks(DVtickLocs, DVtickLabels)
plt.colorbar(label = 'FT Selectivity Index', location = 'bottom', pad = 0.2, shrink = 0.8)
plt.title('FT selectivity map', fontsize = fontSizeTitles)
axColorMapFT.spines["right"].set_visible(False)
axColorMapFT.spines["top"].set_visible(False)


plt.sca(axFtDV)
plt.scatter(bestSelectivityIndexFt[speechResponsive & ~excludeCells], y_coords[speechResponsive & ~excludeCells], c = colorPts,  alpha = 0.5, s = 5)
plt.boxplot(quantilesFT_DV, positions = binsDV, vert = False, widths = 6, boxprops = dict(linewidth = 2), medianprops = dict(linewidth = 2, color = 'k'), showfliers = False)
plt.plot([avgShuffledSI, avgShuffledSI] ,[215,60], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.ylim(215, 60)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('FT Selectivity Index', fontsize = fontSizeLabels)
plt.ylabel('Dorsal - Ventral (mm)', fontsize = fontSizeLabels)
axFtDV.spines["right"].set_visible(False)
axFtDV.spines["top"].set_visible(False)



plt.sca(axFtAP)
plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells], bestSelectivityIndexFt[speechResponsive & ~excludeCells], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesFT_AP, positions = binsAP, vert = True, widths = 2.5, boxprops = dict(linewidth = 2), medianprops = dict(linewidth = 2, color = 'k'), showfliers = False)
plt.plot([165,225], [avgShuffledSI, avgShuffledSI], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.xlabel('Anterior - Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('FT Selectivity Index', fontsize = fontSizeLabels)
axFtAP.spines["right"].set_visible(False)
axFtAP.spines["top"].set_visible(False)

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
