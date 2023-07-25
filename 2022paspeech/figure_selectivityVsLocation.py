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
SAVE_FIGURE = 0
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_selectivityIndices' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 10.5] # In inches
STATSUMMARY = 1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
fontSizePanel = figparams.fontSizePanel
#colorMap = cm.get_cmap('Greens')
#newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
colorMapVOT = cm.get_cmap('Reds')
newMapVOT = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMapVOT(np.linspace(0.2,1,50)))

colorMapFT = cm.get_cmap('Purples')
newMapFT = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMapFT(np.linspace(0.2,1,50)))

colorPts = cp.TangoPalette['SkyBlue1']
colorRandSI = cp.TangoPalette['Aluminium4']
colorVotSelective = cp.TangoPalette['ScarletRed2']
colorFtSelective = cp.TangoPalette['Plum3']
colorNotSelective = cp.TangoPalette['Aluminium3']


figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
z_coords_jittered = figData['z_coords_jittered']
x_coords_jittered = figData['x_coords_jittered']
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
bestSelectivityIndexVot = figData['bestSelectivityIndexVot']
bestSelectivityIndexFt = figData['bestSelectivityIndexFt']
pvalPermutationtestFt = figData['pvalPermutationtestFt']
pvalPermutationtestVot = figData['pvalPermutationtestVot']
shuffledVotBest = figData['shuffledVotBest']
shuffledFtBest = figData['shuffledFtBest']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
avgShuffledSIVot = np.mean(shuffledVotBest[speechResponsive],1)
avgShuffledSIFt = np.mean(shuffledFtBest[speechResponsive], 1)
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



#plt.figure()
gsMain = gridspec.GridSpec(2, 1)

gsVOT = gsMain[0].subgridspec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axColorMapVOT = plt.subplot(gsVOT[:,0])
axDVAPmapVot = gsVOT[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45])
axQuadsVot = plt.subplot(axDVAPmapVot[0,0])
axByAnimalVot = plt.subplot(axDVAPmapVot[1,0])
axVotDonuts = gsVOT[0:,1].subgridspec(4, 1)
axVotAudP = plt.subplot(axVotDonuts[1,0])
axVotAudD = plt.subplot(axVotDonuts[0,0])
axVotAudV = plt.subplot(axVotDonuts[2,0])
axVotTeA = plt.subplot(axVotDonuts[3,0])

gsFT = gsMain[1].subgridspec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axColorMapFT = plt.subplot(gsFT[:,0])
axDVAPmapFt = gsFT[:,2].subgridspec(2,1, height_ratios = [0.55,0.45])
axQuadsFt = plt.subplot(axDVAPmapFt[0,0])
axByAnimalFt = plt.subplot(axDVAPmapFt[1,0])
axFtDonuts = gsFT[0:,1].subgridspec(4, 1)
axFtAudP = plt.subplot(axFtDonuts[1,0])
axFtAudD = plt.subplot(axFtDonuts[0,0])
axFtAudV = plt.subplot(axFtDonuts[2,0])
axFtTeA = plt.subplot(axFtDonuts[3,0])

gsMain.update(left=0.08, right=0.96, top=0.92, bottom=0.08, wspace=0.25, hspace=0.3)
plt.subplots_adjust(top = 0.9, bottom = 0.05, hspace = 0.45, left = 0.05)

nBins = 2
nCompar = np.sum(np.arange(nBins-1,0,-1))

spreadDV = np.max(y_coords[soundResponsive]) - np.min(y_coords[soundResponsive])
binSizeDV = spreadDV/nBins
binsDV = np.arange(np.min(y_coords[soundResponsive]), np.max(y_coords[soundResponsive]), binSizeDV)
binsDV_AAtransform = np.round((binsDV-10)*0.025,1)


spreadAP = np.max(z_coords[soundResponsive]) - np.min(z_coords[soundResponsive])
binSizeAP = spreadAP/nBins
binsAP = np.arange(np.min(z_coords[soundResponsive]), np.max(z_coords[soundResponsive]), binSizeAP)
binsAP_AAtransform = np.round(-0.94 - (280-binsAP)*0.025,1)


votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive
singleSelective = np.logical_xor(votSelective, ftSelective)
mixedSelective = votSelective & ftSelective


# -- QUADRANTS
'''
quantilesDV = np.zeros([nBins, len(soundResponsive)], dtype = bool)
quantilesAP = np.zeros([nBins, len(soundResponsive)], dtype = bool)


audCtxAPbounds = np.array([-3.95,-1.9])
spreadAP = audCtxAPbounds[1]-audCtxAPbounds[0]
quadrantAPThreshold = audCtxAPbounds[0] + (spreadAP/2)
quadrantBoundsAP = np.array([audCtxAPbounds[0], quadrantAPThreshold, audCtxAPbounds[1]])
quadrantBoundsAP_AAtransform = 280 + (quadrantBoundsAP + 0.94)/0.025

audCtxDVbounds = np.array([1.65, 4.65])
spreadDV = audCtxDVbounds[1] - audCtxDVbounds[0]
quadrantDVThreshold = audCtxDVbounds[0] + (spreadDV/2)
quadrantBoundsDV = np.array([audCtxDVbounds[0], quadrantDVThreshold, audCtxDVbounds[1]])
quadrantBoundsDV_AAtransform = (quadrantBoundsDV/0.025)+10


quantilesDV = np.zeros([nBins, len(soundResponsive)], dtype = bool)
quantilesAP = np.zeros([nBins, len(soundResponsive)], dtype = bool)
for indBin, thisBin in enumerate(quantilesDV):
    thisQuantileDV = (y_coords >= quadrantBoundsDV_AAtransform[indBin]) & (y_coords < quadrantBoundsDV_AAtransform[indBin+1])
    thisQuantileAP = (z_coords >= quadrantBoundsAP_AAtransform[indBin]) & (z_coords < quadrantBoundsAP_AAtransform[indBin+1])
    quantilesAP[indBin] = thisQuantileAP
    quantilesDV[indBin] = thisQuantileDV

quadrantLabelsAP = ['posterior', 'anterior']
quadrantLabelsDV = ['dorsal', 'ventral']
quadrantLabels = ['dorsal posterior', 'dorsal anterior', 'ventral posterior', 'ventral anterior']

quadrantTotals = np.array(np.zeros(len(quadrantLabels), dtype=int))
quadrantsVOT = []
quadrantsFT = []
quadrantsVotSelective = []
quadrantsFtSelective = []
quadrantsSingleSelective = []
quadrantsMixedSelective = []
quadrantsSpeechResponsive = []
quadrantsSoundResponsive = []

a = -1
for indBinDV, thisQuantileDV in enumerate(quantilesDV):
    for indBinAP, thisQuantileAP in enumerate(quantilesAP):
        a = a+1
        quadrantTotals[a] = np.sum(thisQuantileDV & thisQuantileAP)
        quadrantsVOT.append(bestSelectivityIndexVot[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsFT.append(bestSelectivityIndexFt[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsVotSelective.append(votSelective[thisQuantileAP & thisQuantileDV & speechResponsive])
        quadrantsFtSelective.append(ftSelective[thisQuantileAP & thisQuantileDV & speechResponsive])
        quadrantsSingleSelective.append(singleSelective[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsMixedSelective.append(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsSpeechResponsive.append(speechResponsive[thisQuantileDV & thisQuantileAP])
        quadrantsSoundResponsive.append(soundResponsive[thisQuantileDV & thisQuantileAP])



# -- TEST EACH ANIMAL
CHECKBYANIMAL = 0
if CHECKBYANIMAL:
    for indMouse, thisMouse in enumerate(np.unique(celldb.subject)):
        print(thisMouse)
        for indBinDV, thisQuantileDV in enumerate(quantilesDV):
            for indBinAP, thisQuantileAP in enumerate(quantilesAP):
                print(f'Quadrant: {quadrantLabelsDV[indBinDV]} {quadrantLabelsAP[indBinAP]}')
                print(f'n Total {np.sum(thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse))}')
                print(f'n SoundResponsive {np.sum(soundResponsive[thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse)])}')
                print(f'n SpeechResponsive {np.sum(speechResponsive[thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse)])}')
                print(f'n Speech Selective {np.sum(singleSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)]) + np.sum(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])}')
                print(f'n Vot Selective {np.sum(votSelective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])}')
                print(f'n Ft Selective {np.sum(ftSelective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])}')
                print(f'n MixedSelective {np.sum(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])}')
'''
#baselineRateByQuad = np.array([celldb.amFiringRateBaseline[quantilesDV[0] & quantilesAP[0]], celldb.amFiringRateBaseline[quantilesDV[0] & quantilesAP[1]], celldb.amFiringRateBaseline[quantilesDV[1] & quantilesAP[0]], celldb.amFiringRateBaseline[quantilesDV[1] & quantilesAP[1]]], dtype=object)


#evokedRateByQuad = np.array([celldb.speechFiringRateBestSustain[quantilesDV[0] & quantilesAP[0]], celldb.speechFiringRateBestSustain[quantilesDV[0] & quantilesAP[1]], celldb.speechFiringRateBestSustain[quantilesDV[1] & quantilesAP[0]], celldb.speechFiringRateBestSustain[quantilesDV[1] & quantilesAP[1]]], dtype = object)



if STATSUMMARY:
    # -- QUADRANTS
    kstat, pvalQuadrantsVOT = stats.kruskal(quadrantsVOT[0], quadrantsVOT[1], quadrantsVOT[2], quadrantsVOT[3])
    kstat, pvalQuadrantsFT = stats.kruskal(quadrantsFT[0], quadrantsFT[1], quadrantsFT[2], quadrantsFT[3])
    nQuadCompar = 6
    quadAlpha = 0.05/6
    quadComparLabels = np.array(['DP vs DA', 'DP vs. VP', 'DP vs. VA', 'DA vs. VP', 'DA vs. VA', 'VP vs. VA'])

    if pvalQuadrantsVOT < 0.05:
        quadrantComparVot = np.ones(nQuadCompar)
        a = -1
        for indBin, thisBin in enumerate(quadrantsVOT):
            nBinCompar = 4 - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quadrantsVOT[indBin], quadrantsVOT[x + indBin + 1])
                quadrantComparVot[a] = pvalMannU
    if pvalQuadrantsFT < 0.05:
        quadrantComparFt = np.ones(nBinCompar)
        a = -1
        for indBin, thisBin in enumerate(quadrantsFT):
            nBinCompar = 4 - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quadrantsFT[indBin], quadrantsFT[x + indBin + 1])
                quadrantComparFt[a] = pvalMannU

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
            oddsratio, pvalFracSpeechResponsive_allcells = stats.fisher_exact(np.array([[np.sum(quadrantsSpeechResponsive[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1])],[len(quadrantsSoundResponsive[indBin]) -np.sum(quadrantsSpeechResponsive[indBin]), len(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsSpeechResponsive[x + indBin + 1])]]))
            quadrantComparFracSpeechResponsive_allcells[a] = pvalFracSpeechResponsive_allcells

    ##--Test frac speech responsive out of soundResponsive
    quadrantComparFracSpeechResponsive_soundResp = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsSpeechResponsive):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracSpeechResponsive_soundResp = stats.fisher_exact(np.array([[np.sum(quadrantsSpeechResponsive[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1])],[np.sum(quadrantsSoundResponsive[indBin]) -np.sum(quadrantsSpeechResponsive[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsSpeechResponsive[x + indBin + 1])]]))
            quadrantComparFracSpeechResponsive_soundResp[a] = pvalFracSpeechResponsive_soundResp

    ##--Test frac VOT selective
    quadrantComparFracVotSelective = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsVotSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsVotSelective[x + indBin + 1])],[np.sum(quadrantsSpeechResponsive[indBin]) -np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1]) - np.sum(quadrantsVotSelective[x + indBin + 1])]]))
            quadrantComparFracVotSelective[a] = pvalFracVOTSelective

    ##--Test frac FT selective
    quadrantComparFracFtSelective = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsFtSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsFtSelective[x + indBin + 1])],[np.sum(quadrantsSpeechResponsive[indBin]) - np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1]) - np.sum(quadrantsFtSelective[x + indBin + 1])]]))
            quadrantComparFracFtSelective[a] = pvalFracFTSelective

    ##--Test frac mixed selective
    quadrantComparFracMixedSelective = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsMixedSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quadrantsSingleSelective[indBin]), np.sum(quadrantsSingleSelective[x + indBin + 1])],[np.sum(quadrantsMixedSelective[indBin]), np.sum(quadrantsMixedSelective[x + indBin + 1])]]))
            quadrantComparFracMixedSelective[a] = pvalFracMixedSelective
    #
    print('--Quadrants Stat Summary --')
    print(f'DP Quadrant: total n = {len(quadrantsSoundResponsive[0])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[0])} ({np.round((np.sum(quadrantsSoundResponsive[0])/len(quadrantsSoundResponsive[0]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[0])} ({np.round((np.sum(quadrantsSpeechResponsive[0])/np.sum(quadrantsSoundResponsive[0]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[0])/len(quadrantsSoundResponsive[0]))*100,1)}% of total), n VOT selective = {np.sum(quadrantsVotSelective[0])} ({np.round((np.sum(quadrantsVotSelective[0])/np.sum(quadrantsSpeechResponsive[0]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[0])} ({np.round((np.sum(quadrantsFtSelective[0])/np.sum(quadrantsSpeechResponsive[0]))*100,1)}% of speech responsive), n mixed selective = {np.sum(quadrantsMixedSelective[0])} ({np.round((np.sum(quadrantsMixedSelective[0])/(np.sum(quadrantsMixedSelective[0]) + np.sum(quadrantsSingleSelective[0])))*100,1)}% of speech selective)')
    print(f'DA Quadrant: total n = {len(quadrantsSoundResponsive[1])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[1])} ({np.round((np.sum(quadrantsSoundResponsive[1])/len(quadrantsSoundResponsive[1]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[1])} ({np.round((np.sum(quadrantsSpeechResponsive[1])/np.sum(quadrantsSoundResponsive[1]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[1])/len(quadrantsSoundResponsive[1]))*100,1)}% of total), n VOT selective = {np.sum(quadrantsVotSelective[1])} ({np.round((np.sum(quadrantsVotSelective[1])/np.sum(quadrantsSpeechResponsive[1]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[1])} ({np.round((np.sum(quadrantsFtSelective[1])/np.sum(quadrantsSpeechResponsive[1]))*100,1)}% of speech responsive), n mixed selective = {np.sum(quadrantsMixedSelective[1])} ({np.round((np.sum(quadrantsMixedSelective[1])/(np.sum(quadrantsMixedSelective[1]) + np.sum(quadrantsSingleSelective[1])))*100,1)}% of speech selective)')
    print(f'VP Quadrant: total n = {len(quadrantsSoundResponsive[2])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[2])} ({np.round((np.sum(quadrantsSoundResponsive[2])/len(quadrantsSoundResponsive[2]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[2])} ({np.round((np.sum(quadrantsSpeechResponsive[2])/np.sum(quadrantsSoundResponsive[2]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[2])/len(quadrantsSoundResponsive[2]))*100,1)}% of total), n VOT selective = {np.sum(quadrantsVotSelective[2])} ({np.round((np.sum(quadrantsVotSelective[2])/np.sum(quadrantsSpeechResponsive[2]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[2])} ({np.round((np.sum(quadrantsFtSelective[2])/np.sum(quadrantsSpeechResponsive[2]))*100,1)}% of speech responsive), n mixed selective = {np.sum(quadrantsMixedSelective[2])} ({np.round((np.sum(quadrantsMixedSelective[2])/(np.sum(quadrantsMixedSelective[2]) + np.sum(quadrantsSingleSelective[2])))*100,1)}% of speech selective)')
    print(f'VA Quadrant: total n = {len(quadrantsSoundResponsive[3])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[3])} ({np.round((np.sum(quadrantsSoundResponsive[3])/len(quadrantsSoundResponsive[3]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[3])} ({np.round((np.sum(quadrantsSpeechResponsive[3])/np.sum(quadrantsSoundResponsive[3]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[3])/len(quadrantsSoundResponsive[3]))*100,1)}% of total), n VOT selective = {np.sum(quadrantsVotSelective[3])} ({np.round((np.sum(quadrantsVotSelective[3])/np.sum(quadrantsSpeechResponsive[3]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[3])} ({np.round((np.sum(quadrantsFtSelective[3])/np.sum(quadrantsSpeechResponsive[3]))*100,1)}% of speech responsive), n mixed selective = {np.sum(quadrantsMixedSelective[3])} ({np.round((np.sum(quadrantsMixedSelective[3])/(np.sum(quadrantsMixedSelective[3]) + np.sum(quadrantsSingleSelective[3])))*100,1)}% of speech selective)')

    if any(quadrantComparFracSoundResponsive < quadAlpha):
        print(f'Significant Quadrant comparisons Frac Sound Responsive: {quadComparLabels[quadrantComparFracSoundResponsive < quadAlpha]}')
    if any(quadrantComparFracSpeechResponsive_allcells < quadAlpha):
        print(f'Significant Quadrant comparisons Frac Speech Responsive out of total cells: {quadComparLabels[quadrantComparFracSpeechResponsive_allcells < quadAlpha]}')
    if any(quadrantComparFracSpeechResponsive_soundResp < quadAlpha):
        print(f'Significant Quadrant comparisons Frac Speech Responsive out of sound responsive cells: {quadComparLabels[quadrantComparFracSpeechResponsive_soundResp < quadAlpha]}')
    if any(quadrantComparFracVotSelective < quadAlpha):
        print(f'Significant Quadrant comparisons Frac VOT Selective: {quadComparLabels[quadrantComparFracVotSelective < quadAlpha]}')
    if any(quadrantComparFracFtSelective < quadAlpha):
        print(f'Significant Quadrant comparisons Frac FT Selective: {quadComparLabels[quadrantComparFracFtSelective < quadAlpha]}')
    if any(quadrantComparFracMixedSelective < quadAlpha):
        print(f'Significant Quadrant comparisons Frac Mixed Selective: {quadComparLabels[quadrantComparFracMixedSelective < quadAlpha]}')



#APtickLocs = np.array([176, 196, 216])
APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
#DVtickLocs = np.array([190, 140, 90])
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
#plt.suptitle('VOT selectivity by location', fontsize = fontSizeTitles)

plt.sca(axColorMapVOT)
#plt.scatter(z_coords_jittered[speechResponsive], y_coords[speechResponsive], c = bestSelectivityIndexVot[speechResponsive], cmap = newMap, s = 3)
plt.scatter(z_coords_jittered[speechResponsive], y_coords[speechResponsive], c = bestSelectivityIndexVot[speechResponsive], cmap = newMapVOT, s = 3)
#plt.ylim(215,60)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
#plt.xlim(165,225)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.colorbar(label = 'VOT Selectivity Index', shrink = 0.8, ticks = [0.2, 0.4, 0.6, 0.8, 1.0])
plt.title('VOT selectivity', fontsize = fontSizeTitles)
axColorMapVOT.spines["right"].set_visible(False)
axColorMapVOT.spines["top"].set_visible(False)

plt.sca(axVotAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0] ))
nSoundResponsiveAudP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[0]))
nTotalAudP = np.sum(recordingAreaName == audCtxAreas[0])
nVOTselectiveAudP = np.sum(votSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nTotalAudP - nVOTselectiveAudP,  nVOTselectiveAudP], colors = [colorNotSelective, colorVotSelective])
axVotAudP.add_artist(circle1)
plt.title(f'AudP,\n n = {int(nTotalAudP)}')

plt.sca(axVotAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1] ))
nSoundResponsiveAudD = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[1]))
nTotalAudD = np.sum(recordingAreaName == audCtxAreas[1])
nVOTselectiveAudD = np.sum(votSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nTotalAudD - nVOTselectiveAudD,  nVOTselectiveAudD], colors = [colorNotSelective, colorVotSelective])
axVotAudD.add_artist(circle2)
plt.title(f'AudD,\n n = {int(nTotalAudD)}')
#plt.legend(labels = ['VOT non-selective', 'VOT selective'], loc = 'upper right', ncol = 2, bbox_to_anchor = (3, 2))

plt.sca(axVotAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2] ))
nSoundResponsiveAudV = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[2]))
nTotalAudV = np.sum(recordingAreaName == audCtxAreas[2])
nVOTselectiveAudV = np.sum(votSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nTotalAudV - nVOTselectiveAudV,  nVOTselectiveAudV], colors = [colorNotSelective, colorVotSelective])
axVotAudV.add_artist(circle3)
plt.title(f'AudV,\n n = {int(nTotalAudV)}')

plt.sca(axVotTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3] ))
nSoundResponsiveTeA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[3] ))
nTotalTeA = np.sum(recordingAreaName == audCtxAreas[3])
nVOTselectiveTeA = np.sum(votSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nTotalTeA - nVOTselectiveTeA,  nVOTselectiveTeA], colors = [colorNotSelective, colorVotSelective])
axVotTeA.add_artist(circle4)
plt.title(f'TeA,\n n = {int(nTotalTeA)}')
plt.legend(labels = ['VOT non-selective', 'VOT selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))

plt.sca(axQuadsVot)
nVOTselectiveDP = np.sum(quadrantsVotSelective[0])
nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0])
nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0])
nTotalDP = quadrantTotals[0]
fracVOTselectiveDP = nVOTselectiveDP/nTotalDP
DPalphaVot = 1

nVOTselectiveDA = np.sum(quadrantsVotSelective[1])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1])
nTotalDA = quadrantTotals[1]
fracVOTselectiveDA = nVOTselectiveDA/nTotalDA
DAalphaVot = fracVOTselectiveDA/fracVOTselectiveDP

nVOTselectiveVP = np.sum(quadrantsVotSelective[2])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2])
nTotalVP = quadrantTotals[2]
fracVOTselectiveVP = nVOTselectiveVP/nTotalVP
VPalphaVot = fracVOTselectiveVP/fracVOTselectiveDP

nVOTselectiveVA = np.sum(quadrantsVotSelective[3])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3])
nTotalVA = quadrantTotals[3]
fracVOTselectiveVA = nVOTselectiveVA/nTotalVA
VAalphaVot = fracVOTselectiveVA/fracVOTselectiveDP


pltVotDP = plt.Rectangle((0,0.5), 0.5, 0.5, facecolor = colorVotSelective, alpha = DPalphaVot)
pltVotDA = plt.Rectangle((0.5,0.5), 0.5, 0.5, facecolor = colorVotSelective, alpha = DAalphaVot)
pltVotVP = plt.Rectangle((0,0), 0.5, 0.5, facecolor = colorVotSelective, alpha = VPalphaVot)
pltVotVA = plt.Rectangle((0.5,0), 0.5, 0.5, facecolor = colorVotSelective, alpha = VAalphaVot)
axQuadsVot.add_artist(pltVotDP)
axQuadsVot.add_artist(pltVotDA)
axQuadsVot.add_artist(pltVotVP)
axQuadsVot.add_artist(pltVotVA)
plt.annotate(f'{np.round(fracVOTselectiveDP*100,1)}% DP\n n = {nTotalDP}', (0.1, 0.65), fontsize = fontSizeTicks, color = 'w')
plt.annotate(f'{np.round(fracVOTselectiveDA*100,1)}% DA\n n = {nTotalDA}', (0.6, 0.65), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracVOTselectiveVP*100,1)}% VP\n n = {nTotalVP}', (0.1, 0.15), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracVOTselectiveVA*100,1)}% VA\n n = {nTotalVA}', (0.6, 0.15), fontsize = fontSizeTicks, color = 'k')
plt.axis('off')




plt.sca(axByAnimalVot)
plt.bar([1, 2, 3, 4], [fracVOTselectiveDP, fracVOTselectiveDA, fracVOTselectiveVA, fracVOTselectiveVP], facecolor = colorVotSelective)
plt.xticks([1,2,3,4], ['DP', 'DA', 'VA', 'VP'])
plt.ylabel('Fraction VOT selective')




plt.sca(axColorMapFT)
plt.scatter(z_coords_jittered[speechResponsive], y_coords[speechResponsive], c = bestSelectivityIndexFt[speechResponsive], cmap = newMapFT, s = 3)
#plt.ylim(215,60)
#plt.xlim(165,225)
plt.ylim(220,40)
plt.xlim(146, 246)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.xticks(APtickLocs, APtickLabels)
plt.yticks(DVtickLocs, DVtickLabels)
plt.colorbar(label = 'FT Selectivity Index', shrink = 0.8, ticks = [0.2, 0.4, 0.6, 0.8, 1.0], extend = 'max', extendrect = True, extendfrac = 0.22)
plt.title('FT selectivity', fontsize = fontSizeTitles)
axColorMapFT.spines["right"].set_visible(False)
axColorMapFT.spines["top"].set_visible(False)



plt.sca(axFtAudP)
circle5 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudP = np.sum(ftSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nTotalAudP - nFTselectiveAudP,  nFTselectiveAudP], colors = [colorNotSelective, colorFtSelective])
axFtAudP.add_artist(circle5)
plt.title(f'AudP,\n n = {int(nTotalAudP)}')

plt.sca(axFtAudD)
circle6 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudD = np.sum(ftSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nTotalAudD - nFTselectiveAudD,  nFTselectiveAudD], colors = [colorNotSelective, colorFtSelective])
axFtAudD.add_artist(circle6)
plt.title(f'AudD,\n n = {int(nTotalAudD)}')
#plt.legend(labels = ['FT non-selective', 'FT selective'], loc = 'upper right', ncol = 2, bbox_to_anchor = (3, 2))


plt.sca(axFtAudV)
circle7 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudV = np.sum(ftSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nTotalAudV - nFTselectiveAudV,  nFTselectiveAudV], colors = [colorNotSelective, colorFtSelective])
axFtAudV.add_artist(circle7)
plt.title(f'AudV,\n n = {int(nTotalAudV)}')

plt.sca(axFtTeA)
circle8 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveTeA = np.sum(ftSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nTotalTeA - nFTselectiveTeA,  nFTselectiveTeA], colors = [colorNotSelective, colorFtSelective])
axFtTeA.add_artist(circle8)
plt.title(f'TeA,\n n = {int(nTotalTeA)}')
plt.legend(labels = ['FT non-selective', 'FT selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))


plt.sca(axQuadsFt)
nFTselectiveDP = np.sum(quadrantsFtSelective[0])
nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0])
nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0])
nTotalDP = quadrantTotals[0]
fracFTselectiveDP = nFTselectiveDP/nTotalDP
DPalphaFt = 1

nFTselectiveDA = np.sum(quadrantsFtSelective[1])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1])
nTotalDA = quadrantTotals[1]
fracFTselectiveDA = nFTselectiveDA/nTotalDA
DAalphaFt = fracFTselectiveDA/fracFTselectiveDP

nFTselectiveVP = np.sum(quadrantsFtSelective[2])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2])
nTotalVP = quadrantTotals[2]
fracFTselectiveVP = nFTselectiveVP/nTotalVP
VPalphaFt = fracFTselectiveVP/fracFTselectiveDP

nFTselectiveVA = np.sum(quadrantsFtSelective[3])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3])
nTotalVA = quadrantTotals[3]
fracFTselectiveVA = nFTselectiveVA/nTotalVA
VAalphaFt = fracFTselectiveVA/fracFTselectiveDP


pltFtDP = plt.Rectangle((0,0.5), 0.5, 0.5, facecolor = colorFtSelective, alpha = DPalphaFt)
pltFtDA = plt.Rectangle((0.5,0.5), 0.5, 0.5, facecolor = colorFtSelective, alpha = DAalphaFt)
pltFtVP = plt.Rectangle((0,0), 0.5, 0.5, facecolor = colorFtSelective, alpha = VPalphaFt)
pltFtVA = plt.Rectangle((0.5,0), 0.5, 0.5, facecolor = colorFtSelective, alpha = VAalphaFt)
axQuadsFt.add_artist(pltFtDP)
axQuadsFt.add_artist(pltFtDA)
axQuadsFt.add_artist(pltFtVP)
axQuadsFt.add_artist(pltFtVA)
plt.annotate(f'{np.round(fracFTselectiveDP*100,1)}% DP\n n = {nTotalDP}', (0.1, 0.65), fontsize = fontSizeTicks, color = 'w')
plt.annotate(f'{np.round(fracFTselectiveDA*100,1)}% DA\n n = {nTotalDA}', (0.6, 0.65), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracFTselectiveVP*100,1)}% VP\n n = {nTotalVP}', (0.1, 0.15), fontsize = fontSizeTicks, color = 'k')
plt.annotate(f'{np.round(fracFTselectiveVA*100,1)}% VA\n n = {nTotalVA}', (0.6, 0.15), fontsize = fontSizeTicks, color = 'k')

plt.axis('off')


plt.sca(axByAnimalFt)
plt.bar([1, 2, 3, 4], [fracFTselectiveDP, fracFTselectiveDA, fracFTselectiveVA, fracFTselectiveVP], facecolor = colorFtSelective)
plt.xticks([1,2,3,4], ['DP', 'DA', 'VA', 'VP'])
plt.ylabel('Fraction FT selective')



labelPosX = [0.05, 0.43, 0.63] # Horiz position for panel labels
labelPosY = [0.94, 0.72,  0.46, 0.24]    # Vert position for panel labels

axColorMapVOT.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axVotAudD.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axQuadsVot.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axByAnimalVot.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axColorMapFT.annotate('E', xy=(labelPosX[0],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axFtAudD.annotate('F', xy=(labelPosX[1],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axQuadsFt.annotate('G', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axByAnimalFt.annotate('H', xy=(labelPosX[2],labelPosY[3]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
