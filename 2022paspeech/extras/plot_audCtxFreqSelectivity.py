import sys
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import celldatabase
from jaratoolbox import colorpalette as cp
sys.path.append('..')
import figparams
import studyparams
import studyutils
from scipy import stats
from importlib import reload

reload(figparams)

FIGNAME = 'selectivityIndices'
SAVE_FIGURE = 1
STATSUMMARY = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'audCtxFreqSelectivity' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [6, 7] # In inches

figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)
databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(databaseDir, 'fulldb_speech_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)

boundDataFile = 'brain_areas_boundaries.npz'
boundDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
boundDataFullPath = os.path.join(boundDataDir,boundDataFile)
boundData = np.load(boundDataFullPath, allow_pickle = True)
contours = boundData['contours']
extentAP = boundData['extentAP']
extentDV = boundData['extentDV']
colorBounds = cp.TangoPalette['Aluminium3'] #'0.75'

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
fontSizeTitles = figparams.fontSizeTitles

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
x_coords_jittered = figData['x_coords_jittered']
z_coords_jittered = figData['z_coords_jittered']
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
toneResponsive = figData['toneResponsive']
bestSelectivityIndexVot = figData['bestSelectivityIndexVot']
bestSelectivityIndexFt = figData['bestSelectivityIndexFt']
isCortical = figData['isCortical']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
quantilesDV = figData['quantilesDV']
quantilesAP = figData['quantilesAP']
quadrantLabels = figData['quadrantLabels']
quadrantsVOT = figData['quadrantsVOT']
quadrantsFT = figData['quadrantsFT']
quadrantsVOT = figData['quadrantsVOT']
quadrantsFT = figData['quadrantsFT']
quadrantsVotSelective = figData['quadrantsVotSelective']
quadrantsFtSelective = figData['quadrantsFtSelective']
pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
possibleFreqs= np.unique(celldb.toneCharactFreq[~np.isnan(celldb.toneCharactFreq)])
votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive


APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)

logFreq = np.log2(celldb.toneCharactFreq[isCortical])
possibleLogFreq = np.log2(possibleFreqs)
possibleLogFreqStrings = [f'{(2**freq)/1000:0.1f}' for freq in possibleLogFreq]

plt.clf()
ax = plt.gca()
#cmap = 'viridis'
cmap = 'brg_r'
#cmap = 'rainbow'
plt.scatter(z_coords_jittered[isCortical], y_coords[isCortical], c=logFreq, cmap=cmap, s=8)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
#cbar = plt.colorbar()
cbar = plt.colorbar(ticks = possibleLogFreq)
cbar.set_ticklabels(possibleLogFreqStrings)
cbar.set_label('Best Frequency (kHz)', rotation=270, labelpad=12)
#cbar.ax.set_yticklabels(np.round(possibleFreqs/1000,1))
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_aspect('equal')
studyutils.plot_quadrants(ax, extentAP, extentDV, color=colorBounds)
for contour in contours:
    plt.plot(contour[:, 1], contour[:, 0], linewidth=1.5, color=colorBounds, clip_on=False, zorder=-1)
plt.tight_layout()
plt.show()

extraplots.save_figure(figFilename, figFormat, figSize, outputDir)


if STATSUMMARY:
    nQuadCompar = 6
    atlasAreasBestFreq = []
    for indArea, thisArea in enumerate(audCtxAreas):
        atlasAreasBestFreq.append(celldb.toneCharactFreq[(recordingAreaName == audCtxAreas[indArea]) & speechResponsive])
        #atlasAreasBestFreq.append(celldb.toneCharactFreq[(recordingAreaName == audCtxAreas[indArea]) & toneResponsive])
        #atlasAreasBestFreq[indArea] = atlasAreasBestFreq[indArea][~np.isnan(atlasAreasBestFreq[indArea])]

    quadrantsBestFreq = []
    a = -1
    for indBinDV, thisQuantileDV in enumerate(quantilesDV):
        for indBinAP, thisQuantileAP in enumerate(quantilesAP):
            a = a+1
            quadrantsBestFreq.append(celldb.toneCharactFreq[thisQuantileDV & thisQuantileAP &speechResponsive])
            #quadrantsBestFreq.append(celldb.toneCharactFreq[thisQuantileDV & thisQuantileAP &toneResponsive])
            #quadrantsBestFreq[a] = quadrantsBestFreq[a][~np.isnan(quadrantsBestFreq[a])]

    kstat, pvalAreasBestFreq = stats.kruskal(atlasAreasBestFreq[0][~np.isnan(atlasAreasBestFreq[0])], atlasAreasBestFreq[1][~np.isnan(atlasAreasBestFreq[1])], atlasAreasBestFreq[2][~np.isnan(atlasAreasBestFreq[2])], atlasAreasBestFreq[3][~np.isnan(atlasAreasBestFreq[3])])

    kstat, pvalQuadrantsBestFreq = stats.kruskal(quadrantsBestFreq[0][~np.isnan(quadrantsBestFreq[0])], quadrantsBestFreq[1][~np.isnan(quadrantsBestFreq[1])], quadrantsBestFreq[2][~np.isnan(quadrantsBestFreq[2])], quadrantsBestFreq[3][~np.isnan(quadrantsBestFreq[3])])



    quadComparLabels = np.array(['DP vs DA', 'DP vs. VP', 'DP vs. VA', 'DA vs. VP', 'DA vs. VA', 'VP vs. VA'])
    areaComparLabels = np.array(['AudP vs AudD', 'AudP vs. AudV', 'AudP vs. TeA', 'AudD vs. AudV', 'AudD vs. TeA', 'AudV vs. TeA'])
    if pvalAreasBestFreq < 0.05:
        atlasAreasComparBestFreq = np.ones(nQuadCompar)
        a = -1
        for indBin, thisBin in enumerate(atlasAreasBestFreq):
            nBinCompar = 4 - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(atlasAreasBestFreq[indBin][~np.isnan(atlasAreasBestFreq[indBin])], atlasAreasBestFreq[x + indBin + 1][~np.isnan(atlasAreasBestFreq[x + indBin + 1])])
                atlasAreasComparBestFreq[a] = pvalMannU

    if pvalQuadrantsBestFreq < 0.05:
        quadrantComparBestFreq = np.ones(nQuadCompar)
        a = -1
        for indBin, thisBin in enumerate(quadrantsBestFreq):
            nBinCompar = 4 - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quadrantsBestFreq[indBin][~np.isnan(quadrantsBestFreq[indBin])], quadrantsBestFreq[x + indBin + 1][~np.isnan(quadrantsBestFreq[x + indBin + 1])])
                quadrantComparBestFreq[a] = pvalMannU




    print('--STAT SUMMARY')
    print(f'distribution of best frequencies among areas: p = {pvalAreasBestFreq}')
    print(f'distribution of best frequencies among regions: p = {pvalQuadrantsBestFreq}')
    if pvalAreasBestFreq < 0.05:
        print('atlas-area compare best freqs:')
        for indCompar, thisCompar in enumerate(areaComparLabels):
            print(f'{areaComparLabels[indCompar]} p = {atlasAreasComparBestFreq[indCompar]}')
    if pvalQuadrantsBestFreq < 0.05:
        print('Quadrant compare best freqs:')
        for indCompar, thisCompar in enumerate(quadComparLabels):
            print(f'{quadComparLabels[indCompar]} p = {quadrantComparBestFreq[indCompar]}')
    print('atlas-areas')
    for indArea, thisArea in enumerate(audCtxAreas):
        print(f'{thisArea}: best freq (Hz)  median = {int(np.nanmedian(atlasAreasBestFreq[indArea]))}, mean = {int(np.nanmean(atlasAreasBestFreq[indArea]))}')
    print('quadrant regions')
    for indArea, thisArea in enumerate(quadrantLabels):
        print(f'{thisArea}: best freq (Hz)  median = {int(np.nanmedian(quadrantsBestFreq[indArea]))}, mean = {int(np.nanmean(quadrantsBestFreq[indArea]))}')
