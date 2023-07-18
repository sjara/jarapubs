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
colorMap = cm.get_cmap('Greens')
newMap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n='BuGn', a = 0.1, b = 1), colorMap(np.linspace(0.2,1,50)))
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


#plt.figure()
gsMain = gridspec.GridSpec(2, 1)

gsVOT = gsMain[0].subgridspec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axColorMapVOT = plt.subplot(gsVOT[:,0])
axDVAPmapVot = gsVOT[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45])
axVotDV = plt.subplot(axDVAPmapVot[0,0])
axVotAP = plt.subplot(axDVAPmapVot[1,0])
axVotDonuts = gsVOT[0:,1].subgridspec(4, 1)
axVotAudP = plt.subplot(axVotDonuts[1,0])
axVotAudD = plt.subplot(axVotDonuts[0,0])
axVotAudV = plt.subplot(axVotDonuts[2,0])
axVotTeA = plt.subplot(axVotDonuts[3,0])

gsFT = gsMain[1].subgridspec(4, 3, width_ratios = [0.4, 0.25, 0.35])
axColorMapFT = plt.subplot(gsFT[:,0])
axDVAPmapFt = gsFT[:,2].subgridspec(2,1, height_ratios = [0.55,0.45])
axFtDV = plt.subplot(axDVAPmapFt[0,0])
axFtAP = plt.subplot(axDVAPmapFt[1,0])
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

'''
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
quantilesDV = np.zeros([nBins, len(speechResponsive)], dtype = bool)
quantilesAP = np.zeros([nBins, len(speechResponsive)], dtype = bool)



# -- BINS in AP and DV
for indBin, thisBin in enumerate(binsDV):
    if indBin < len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords < binsDV[indBin+1])
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords < binsAP[indBin+1])
    elif indBin == len(binsDV) - 1:
        thisQuantileDV = (y_coords >= binsDV[indBin]) & (y_coords <= np.max(y_coords[speechResponsive]))
        thisQuantileAP = (z_coords >= binsAP[indBin]) & (z_coords <= np.max(z_coords[speechResponsive]))
    quantilesAP[indBin] = thisQuantileAP
    quantilesDV[indBin] = thisQuantileDV
    quantilesVOT_DV.append(bestSelectivityIndexVot[thisQuantileDV & speechResponsive])
    quantilesFT_DV.append(bestSelectivityIndexFt[thisQuantileDV & speechResponsive])
    quantilesVOT_AP.append(bestSelectivityIndexVot[thisQuantileAP & speechResponsive])
    quantilesFT_AP.append(bestSelectivityIndexFt[thisQuantileAP & speechResponsive])
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
'''

# -- QUADRANTS

quantilesDV = np.zeros([nBins, len(soundResponsive)], dtype = bool)
quantilesAP = np.zeros([nBins, len(soundResponsive)], dtype = bool)
'''
## --OLD--
spreadDV = np.max(y_coords[soundResponsive]) - np.min(y_coords[soundResponsive])
binSizeDV = spreadDV/nBins
binsDV = np.arange(np.min(y_coords[soundResponsive]), np.max(y_coords[soundResponsive]), binSizeDV)
binsDV_AAtransform = np.round((binsDV-10)*0.025,1)


spreadAP = np.max(z_coords[soundResponsive]) - np.min(z_coords[soundResponsive])
binSizeAP = spreadAP/nBins
binsAP = np.arange(np.min(z_coords[soundResponsive]), np.max(z_coords[soundResponsive]), binSizeAP)
binsAP_AAtransform = np.round(-0.94 - (280-binsAP)*0.025,1)
'''

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

quadrantsVOT = []
quadrantsFT = []
quadrantsVotSelective = []
quadrantsFtSelective = []
quadrantsSingleSelective = []
quadrantsMixedSelective = []
quadrantsSpeechResponsive = []
quadrantsSoundResponsive = []


for indBinDV, thisQuantileDV in enumerate(quantilesDV):
    for indBinAP, thisQuantileAP in enumerate(quantilesAP):
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


    '''
    # -- BINS
    binAlpha = 0.05/nCompar
    if nBins == 10:
        kstat, pvalKruskalVotDV= stats.kruskal(quantilesVOT_DV[0], quantilesVOT_DV[1], quantilesVOT_DV[2], quantilesVOT_DV[3], quantilesVOT_DV[4], quantilesVOT_DV[5], quantilesVOT_DV[6], quantilesVOT_DV[7], quantilesVOT_DV[8], quantilesVOT_DV[9])
        kstat, pvalKruskalVotAP= stats.kruskal(quantilesVOT_AP[0], quantilesVOT_AP[1], quantilesVOT_AP[2], quantilesVOT_AP[3], quantilesVOT_AP[4], quantilesVOT_AP[5], quantilesVOT_AP[6], quantilesVOT_AP[7], quantilesVOT_AP[8], quantilesVOT_AP[9])
        kstat, pvalKruskalFtDV= stats.kruskal(quantilesFT_DV[0], quantilesFT_DV[1], quantilesFT_DV[2], quantilesFT_DV[3], quantilesFT_DV[4], quantilesFT_DV[5], quantilesFT_DV[6], quantilesFT_DV[7], quantilesFT_DV[8], quantilesFT_DV[9])
        kstat, pvalKruskalFtAP= stats.kruskal(quantilesFT_AP[0], quantilesFT_AP[1], quantilesFT_AP[2], quantilesFT_AP[3], quantilesFT_AP[4], quantilesFT_AP[5], quantilesFT_AP[6], quantilesFT_AP[7], quantilesFT_AP[8], quantilesFT_AP[9])
    elif nBins == 9:
        kstat, pvalKruskalVotDV= stats.kruskal(quantilesVOT_DV[0], quantilesVOT_DV[1], quantilesVOT_DV[2], quantilesVOT_DV[3], quantilesVOT_DV[4], quantilesVOT_DV[5], quantilesVOT_DV[6], quantilesVOT_DV[7], quantilesVOT_DV[8])
        kstat, pvalKruskalVotAP= stats.kruskal(quantilesVOT_AP[0], quantilesVOT_AP[1], quantilesVOT_AP[2], quantilesVOT_AP[3], quantilesVOT_AP[4], quantilesVOT_AP[5], quantilesVOT_AP[6], quantilesVOT_AP[7], quantilesVOT_AP[8])
        kstat, pvalKruskalFtDV= stats.kruskal(quantilesFT_DV[0], quantilesFT_DV[1], quantilesFT_DV[2], quantilesFT_DV[3], quantilesFT_DV[4], quantilesFT_DV[5], quantilesFT_DV[6], quantilesFT_DV[7], quantilesFT_DV[8])
        kstat, pvalKruskalFtAP= stats.kruskal(quantilesFT_AP[0], quantilesFT_AP[1], quantilesFT_AP[2], quantilesFT_AP[3], quantilesFT_AP[4], quantilesFT_AP[5], quantilesFT_AP[6], quantilesFT_AP[7], quantilesFT_AP[8])
    elif nBins == 8:
        kstat, pvalKruskalVotDV= stats.kruskal(quantilesVOT_DV[0], quantilesVOT_DV[1], quantilesVOT_DV[2], quantilesVOT_DV[3], quantilesVOT_DV[4], quantilesVOT_DV[5], quantilesVOT_DV[6], quantilesVOT_DV[7])
        kstat, pvalKruskalVotAP= stats.kruskal(quantilesVOT_AP[0], quantilesVOT_AP[1], quantilesVOT_AP[2], quantilesVOT_AP[3], quantilesVOT_AP[4], quantilesVOT_AP[5], quantilesVOT_AP[6], quantilesVOT_AP[7])
        kstat, pvalKruskalFtDV= stats.kruskal(quantilesFT_DV[0], quantilesFT_DV[1], quantilesFT_DV[2], quantilesFT_DV[3], quantilesFT_DV[4], quantilesFT_DV[5], quantilesFT_DV[6], quantilesFT_DV[7])
        kstat, pvalKruskalFtAP= stats.kruskal(quantilesFT_AP[0], quantilesFT_AP[1], quantilesFT_AP[2], quantilesFT_AP[3], quantilesFT_AP[4], quantilesFT_AP[5], quantilesFT_AP[6], quantilesFT_AP[7])

    if pvalKruskalVotAP < 0.025:
        binByBinComparVotAP = np.ones(nCompar)
        a = -1
        for indBin, thisBin in enumerate(quantilesVOT_AP):
            nBinCompar = nBins - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quantilesVOT_AP[indBin], quantilesVOT_AP[x + indBin + 1])
                binByBinComparVotAP[a] = pvalMannU

    if pvalKruskalVotDV < 0.025:
        binByBinComparVotDV = np.ones(nCompar)
        a = -1
        for indBin, thisBin in enumerate(quantilesVOT_AP):
            nBinCompar = nBins - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quantilesVOT_DV[indBin], quantilesVOT_DV[x + indBin + 1])
                binByBinComparVotDV[a] = pvalMannU

    if pvalKruskalFtAP < 0.025:
        binByBinComparFtAP = np.ones(nCompar)
        a = -1
        for indBin, thisBin in enumerate(quantilesVOT_AP):
            nBinCompar = nBins - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quantilesFT_AP[indBin], quantilesFT_AP[x + indBin + 1])
                binByBinComparFtAP[a] = pvalMannU

    if pvalKruskalFtDV < 0.025:
        binByBinComparFtDV = np.ones(nCompar)
        a = -1
        for indBin, thisBin in enumerate(quantilesVOT_AP):
            nBinCompar = nBins - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quantilesFT_DV[indBin], quantilesFT_DV[x + indBin + 1])
                binByBinComparFtDV[a] = pvalMannU


    ##--Test frac VOT selective DV
    binByBinComparFracVotSelectiveDV = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quantilesVotSelective_DV[indBin]), np.sum(quantilesVotSelective_DV[x + indBin + 1])],[np.sum(quantilesResponsive_DV[indBin]), np.sum(quantilesResponsive_DV[x + indBin + 1])]]))
            binByBinComparFracVotSelectiveDV[a] = pvalFracVOTSelective

    ##--Test frac VOT selective AP
    binByBinComparFracVotSelectiveAP = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quantilesVotSelective_AP[indBin]), np.sum(quantilesVotSelective_AP[x + indBin + 1])],[np.sum(quantilesResponsive_AP[indBin]), np.sum(quantilesResponsive_AP[x + indBin + 1])]]))
            binByBinComparFracVotSelectiveAP[a] = pvalFracVOTSelective

    ##--Test frac FT selective DV
    binByBinComparFracFtSelectiveDV = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quantilesFtSelective_DV[indBin]), np.sum(quantilesFtSelective_DV[x + indBin + 1])],[np.sum(quantilesResponsive_DV[indBin]), np.sum(quantilesResponsive_DV[x + indBin + 1])]]))
            binByBinComparFracFtSelectiveDV[a] = pvalFracFTSelective


    ##--Test frac FT selective AP
    binByBinComparFracFtSelectiveAP =  np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quantilesFtSelective_AP[indBin]), np.sum(quantilesFtSelective_AP[x + indBin + 1])],[np.sum(    quantilesResponsive_AP[indBin]), np.sum(quantilesResponsive_AP[x + indBin + 1])]]))
            binByBinComparFracFtSelectiveAP[a] = pvalFracFTSelective


    ##--Test frac mixed selective DV
    binByBinComparFracMixedSelectiveDV = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quantilesSingleSelective_DV[indBin]), np.sum(quantilesSingleSelective_DV[x + indBin + 1])],[np.sum(quantilesMixedSelective_DV[indBin]), np.sum(quantilesMixedSelective_DV[x + indBin + 1])]]))
            binByBinComparFracMixedSelectiveDV[a] = pvalFracMixedSelective


    ##--Test frac mixed selective AP
    binByBinComparFracMixedSelectiveAP = np.ones(nCompar)
    a = -1
    for indBin, thisBin in enumerate(quantilesVOT_AP):
        nBinCompar = nBins - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracMixedSelective = stats.fisher_exact(np.array([[np.sum(quantilesSingleSelective_AP[indBin]), np.sum(quantilesSingleSelective_AP[x + indBin + 1])],[np.sum(quantilesMixedSelective_AP[indBin]), np.sum(quantilesMixedSelective_AP[x + indBin + 1])]]))
            binByBinComparFracMixedSelectiveAP[a] = pvalFracMixedSelective

    if nBins == 9:
        binByBinLabels = np.array(['bin0v1', 'bin0v2', 'bin0v3', 'bin0v4', 'bin0v5', 'bin0v6', 'bin0v7', 'bin0v8', 'bin1v2', 'bin1v3', 'bin1v4', 'bin1v5', 'bin1v6', 'bin1v7', 'bin1v8', 'bin2v3', 'bin2v4', 'bin2v5', 'bin2v6', 'bin2v7', 'bin2v8', 'bin3v4', 'bin3v5', 'bin3v6', 'bin3v7', 'bin3v8', 'bin4v5', 'bin4v6', 'bin4v7', 'bin4v8', 'bin5v6', 'bin5v7', 'bin5v8', 'bin6v7', 'bin6v8', 'bin7v8'])
    elif nBins == 10:
        binByBinLabels = np.array(['bin0v1', 'bin0v2', 'bin0v3', 'bin0v4', 'bin0v5', 'bin0v6', 'bin0v7', 'bin0v8', 'bin0v9', 'bin1v2', 'bin1v3', 'bin1v4', 'bin1v5', 'bin1v6', 'bin1v7', 'bin1v8', 'bin1v9', 'bin2v3', 'bin2v4', 'bin2v5', 'bin2v6', 'bin2v7', 'bin2v8', 'bin2v9', 'bin3v4', 'bin3v5', 'bin3v6', 'bin3v7', 'bin3v8', 'bin3v9', 'bin4v5', 'bin4v6', 'bin4v7', 'bin4v8', 'bin4v9', 'bin5v6', 'bin5v7', 'bin5v8', 'bin5v9', 'bin6v7', 'bin6v8', 'bin6v9', 'bin7v8', 'bin7v9', 'bin8v9'])
    elif nBins == 8:
        binByBinLabels = np.array(['bin0v1', 'bin0v2', 'bin0v3', 'bin0v4', 'bin0v5', 'bin0v6', 'bin0v7', 'bin1v2', 'bin1v3', 'bin1v4', 'bin1v5', 'bin1v6', 'bin1v7', 'bin2v3', 'bin2v4', 'bin2v5', 'bin2v6', 'bin2v7', 'bin3v4', 'bin3v5', 'bin3v6', 'bin3v7', 'bin4v5', 'bin4v6', 'bin4v7', 'bin5v6', 'bin5v7', 'bin6v7'])

    print('--Stats Summary--')
    print(f'nBins = {nBins}, nCompar = {nCompar}')
    print(f'pvalKruskal Vot_AP = {np.round(pvalKruskalVotAP, 3)}')
    print(f'pvalKruskal Vot_DV = {np.round(pvalKruskalVotDV, 3)}')
    print(f'pvalKruskal Ft_AP = {np.round(pvalKruskalFtAP, 3)}')
    print(f'pvalKruskal Ft_DV = {np.round(pvalKruskalFtDV,3)}')
    if pvalKruskalVotAP < 0.025:
        print(f'Significant Bin comparisons VOT AP: {binByBinLabels[binByBinComparVotAP < binAlpha]}')
    if pvalKruskalVotDV < 0.025:
        print(f'Significant Bin comparisons VOT DV: {binByBinLabels[binByBinComparVotDV < binAlpha]}')
    if pvalKruskalFtAP < 0.025:
        print(f'Significant Bin comparisons FT AP: {binByBinLabels[binByBinComparFtAP < binAlpha]}')
    if pvalKruskalFtDV < 0.025:
        print(f'Significant Bin comparisons FT DV: {binByBinLabels[binByBinComparFtDV < binAlpha]}')
    if any(binByBinComparFracVotSelectiveAP < binAlpha):
        print(f'Significant Bin comparisons Frac VOT selective AP: {binByBinLabels[binByBinComparFracVotSelectiveAP < binAlpha]}')
    if any(binByBinComparFracVotSelectiveDV < binAlpha):
        print(f'Significant Bin comparisons Frac VOT selective DV: {binByBinLabels[binByBinComparFracVotSelectiveDV < binAlpha]}')
    if any(binByBinComparFracFtSelectiveAP < binAlpha):
        print(f'Significant Bin comparisons Frac FT selective AP: {binByBinLabels[binByBinComparFracFtSelectiveAP < binAlpha]}')
    if any(binByBinComparFracVotSelectiveDV < binAlpha):
        print(f'Significant Bin comparisons Frac FT selective DV: {binByBinLabels[binByBinComparFracFtSelectiveDV < binAlpha]}')
    if any(binByBinComparFracMixedSelectiveAP < binAlpha):
        print(f'Significant Bin comparisons Frac Mixed selective AP: {binByBinLabels[binByBinComparFracMixedSelectiveAP < binAlpha]}')
    if any(binByBinComparFracMixedSelectiveDV < binAlpha):
        print(f'Significant Bin comparisons Frac Mixed selective DV: {binByBinLabels[binByBinComparFracMixedSelectiveDV < binAlpha]}')
    '''
    '''
    pvalFtvsVot_AP = np.ones(len(quantilesFT_AP))
    pvalFtvsVot_DV = np.ones(len(quantilesFT_DV))

    for indBin, thisBin in enumerate(quantilesVOT_AP):
        wstat, pvalFtvsVot_AP[indBin] = stats.wilcoxon(quantilesFT_AP[indBin], quantilesVOT_AP[indBin])
        wstat, pvalFtvsVot_DV[indBin] = stats.wilcoxon(quantilesFT_DV[indBin], quantilesVOT_DV[indBin])
    '''

#APtickLocs = np.array([176, 196, 216])
APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
#DVtickLocs = np.array([190, 140, 90])
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
#plt.suptitle('VOT selectivity by location', fontsize = fontSizeTitles)

plt.sca(axColorMapVOT)
plt.scatter(z_coords_jittered[speechResponsive], y_coords[speechResponsive], c = bestSelectivityIndexVot[speechResponsive], cmap = newMap, s = 3)
#plt.ylim(215,60)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
#plt.xlim(165,225)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.colorbar(label = 'VOT Selectivity Index', shrink = 0.8, ticks = [0.2, 0.4, 0.6, 0.8, 1.0])
plt.title('VOT selectivity', fontsize = fontSizeTitles)
axColorMapVOT.spines["right"].set_visible(False)
axColorMapVOT.spines["top"].set_visible(False)



#plt.annotate('D', xy = (170, -0.75), xycoords = ('data', 'axes fraction'),  xytext = (160 ,-0.55), textcoords = ('data', 'axes fraction'), arrowprops=dict(arrowstyle="<->", connectionstyle = "angle, angleA = -90, angleB = 0, rad=0", lw = 2), fontsize = fontSizeLabels, fontweight = 'bold')
#plt.annotate('A', xy =(170 ,-0.77), xycoords = ('data', 'axes fraction'), fontsize = fontSizeLabels, fontweight = 'bold')
chanceAPposition = [binsAP[-1] + 2.25*binSizeAP]
chanceDVposition = [binsDV[0] - 2*binSizeDV]


'''
plt.sca(axVotDV)
plt.scatter(bestSelectivityIndexVot[speechResponsive], y_coords[speechResponsive], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesVOT_DV, positions = binsDV, vert = False, widths = binSizeDV/1.75, boxprops = dict(linewidth = 1.5), medianprops = dict(linewidth = 1.5, color = 'k'), whiskerprops = dict(linewidth = 1.5, color = 'k'), capprops = dict(linewidth = 1.5, color = 'k'), showfliers = False)
plt.boxplot(avgShuffledSIVot, positions = chanceDVposition, vert = False, widths = binSizeDV/1.75, boxprops = dict(linewidth = 1.5, color = colorRandSI), medianprops = dict(linewidth = 1.5, color = colorRandSI), whiskerprops = dict(linewidth = 1.5, color = colorRandSI), capprops = dict(linewidth = 1.5, color = colorRandSI), showfliers = False )
#plt.plot([avgShuffledSIVot, avgShuffledSIVot] ,[215,60], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.ylim(190, 40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('VOT Selectivity Index', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
axVotDV.spines["right"].set_visible(False)
axVotDV.spines["top"].set_visible(False)



plt.sca(axVotAP)
plt.scatter(z_coords_jittered[speechResponsive], bestSelectivityIndexVot[speechResponsive], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesVOT_AP, positions = binsAP, vert = True, widths = binSizeAP/1.75, boxprops = dict(linewidth = 1.5), medianprops = dict(linewidth = 1.5, color = 'k'), whiskerprops = dict(linewidth = 1.5, color = 'k'), capprops = dict(linewidth = 1.5, color = 'k'), showfliers = False)
plt.boxplot(avgShuffledSIVot, positions = chanceAPposition, vert = True, widths = binSizeAP/1.75, boxprops = dict(linewidth = 1.5, color = colorRandSI), medianprops = dict(linewidth = 1.5, color = colorRandSI), whiskerprops = dict(linewidth = 1.5, color = colorRandSI), capprops = dict(linewidth = 1.5, color = colorRandSI), showfliers = False )
#plt.plot([165,225], [avgShuffledSIVot, avgShuffledSIVot], linestyle = '--', linewidth = 2, c =colorRandSI)
#plt.ylim(215, 60)
#plt.yticks([])
#ax2.invert_yaxis()
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.ylabel('VOT Selectivity Index', fontsize = fontSizeLabels)
axVotAP.spines["right"].set_visible(False)
axVotAP.spines["top"].set_visible(False)
'''


plt.sca(axVotAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0] ))
nVOTselectiveAudP = np.sum(votSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nSpeechResponsiveAudP - nVOTselectiveAudP,  nVOTselectiveAudP], colors = [colorNotSelective, colorVotSelective])
axVotAudP.add_artist(circle1)
plt.title(f'AudP,\n n = {int(nSpeechResponsiveAudP)}')

plt.sca(axVotAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1] ))
nVOTselectiveAudD = np.sum(votSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nSpeechResponsiveAudD - nVOTselectiveAudD,  nVOTselectiveAudD], colors = [colorNotSelective, colorVotSelective])
axVotAudD.add_artist(circle2)
plt.title(f'AudD,\n n = {int(nSpeechResponsiveAudD)}')
#plt.legend(labels = ['VOT non-selective', 'VOT selective'], loc = 'upper right', ncol = 2, bbox_to_anchor = (3, 2))




plt.sca(axVotAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2] ))
nVOTselectiveAudV = np.sum(votSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nSpeechResponsiveAudV - nVOTselectiveAudV,  nVOTselectiveAudV], colors = [colorNotSelective, colorVotSelective])
axVotAudV.add_artist(circle3)
plt.title(f'AudV,\n n = {int(nSpeechResponsiveAudV)}')

plt.sca(axVotTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3] ))
nVOTselectiveTeA = np.sum(votSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nSpeechResponsiveTeA - nVOTselectiveTeA,  nVOTselectiveTeA], colors = [colorNotSelective, colorVotSelective])
axVotTeA.add_artist(circle4)
plt.title(f'TeA,\n n = {int(nSpeechResponsiveTeA)}')
plt.legend(labels = ['VOT non-selective', 'VOT selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))




#plt.suptitle('FT selectivity by location', fontsize = fontSizeTitles)

plt.sca(axColorMapFT)
plt.scatter(z_coords_jittered[speechResponsive], y_coords[speechResponsive], c = bestSelectivityIndexFt[speechResponsive], cmap = newMap, s = 3)
#plt.ylim(215,60)
#plt.xlim(165,225)
plt.ylim(220,40)
plt.xlim(146, 246)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.xticks(APtickLocs, APtickLabels)
plt.yticks(DVtickLocs, DVtickLabels)
plt.colorbar(label = 'FT Selectivity Index', shrink = 0.8, ticks = [0.2, 0.4, 0.6, 0.8, 1.0], extend = 'max', extendrect = True, extendfrac = 0.22)
plt.title('FT selectivity', fontsize = fontSizeTitles)
axColorMapFT.spines["right"].set_visible(False)
axColorMapFT.spines["top"].set_visible(False)

'''
plt.sca(axFtDV)
plt.scatter(bestSelectivityIndexFt[speechResponsive], y_coords[speechResponsive], c = colorPts,  alpha = 0.5, s = 5)
plt.boxplot(quantilesFT_DV, positions = binsDV, vert = False, widths = binSizeDV/1.75, boxprops = dict(linewidth = 1.5), medianprops = dict(linewidth = 1.5, color = 'k'), whiskerprops = dict(linewidth = 1.5, color = 'k'), capprops = dict(linewidth = 1.5, color = 'k'), showfliers = False)
plt.boxplot(avgShuffledSIFt, positions = chanceDVposition, vert = False, widths = binSizeDV/1.75, boxprops = dict(linewidth = 1.5, color = colorRandSI), medianprops = dict(linewidth = 1.5, color = colorRandSI), whiskerprops = dict(linewidth = 1.5, color = colorRandSI), capprops = dict(linewidth = 1.5, color = colorRandSI), showfliers = False )
#plt.plot([avgShuffledSIFt, avgShuffledSIFt] ,[215,60], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.ylim(190, 40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('FT Selectivity Index', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
axFtDV.spines["right"].set_visible(False)
axFtDV.spines["top"].set_visible(False)



plt.sca(axFtAP)
plt.scatter(z_coords_jittered[speechResponsive], bestSelectivityIndexFt[speechResponsive], c = colorPts, alpha = 0.7, s = 5)
plt.boxplot(quantilesFT_AP, positions = binsAP, vert = True, widths = binSizeAP/1.75, boxprops = dict(linewidth = 1.5), medianprops = dict(linewidth = 1.5, color = 'k'), whiskerprops = dict(linewidth = 1.5, color = 'k'), capprops = dict(linewidth = 1.5, color = 'k'), showfliers = False)
plt.boxplot(avgShuffledSIFt, positions = chanceAPposition, vert = True, widths = binSizeAP/1.75, boxprops = dict(linewidth = 1.5, color = colorRandSI, ), medianprops = dict(linewidth = 1.5, color = colorRandSI, ), whiskerprops = dict(linewidth = 1.5, color = colorRandSI, ), capprops = dict(linewidth = 1.5, color = colorRandSI), showfliers = False )
#plt.plot([165,225], [avgShuffledSIFt, avgShuffledSIFt], linestyle = '--', linewidth = 2, c =colorRandSI)
plt.xlim(165,225)
plt.xticks(APtickLocs, APtickLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.ylabel('FT Selectivity Index', fontsize = fontSizeLabels)
axFtAP.spines["right"].set_visible(False)
axFtAP.spines["top"].set_visible(False)
'''

plt.sca(axFtAudP)
circle5 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudP = np.sum(ftSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nSpeechResponsiveAudP - nFTselectiveAudP,  nFTselectiveAudP], colors = [colorNotSelective, colorFtSelective])
axFtAudP.add_artist(circle5)
plt.title(f'AudP,\n n = {int(nSpeechResponsiveAudP)}')

plt.sca(axFtAudD)
circle6 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudD = np.sum(ftSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nSpeechResponsiveAudD - nFTselectiveAudD,  nFTselectiveAudD], colors = [colorNotSelective, colorFtSelective])
axFtAudD.add_artist(circle6)
plt.title(f'AudD,\n n = {int(nSpeechResponsiveAudD)}')
#plt.legend(labels = ['FT non-selective', 'FT selective'], loc = 'upper right', ncol = 2, bbox_to_anchor = (3, 2))


plt.sca(axFtAudV)
circle7 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudV = np.sum(ftSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nSpeechResponsiveAudV - nFTselectiveAudV,  nFTselectiveAudV], colors = [colorNotSelective, colorFtSelective])
axFtAudV.add_artist(circle7)
plt.title(f'AudV,\n n = {int(nSpeechResponsiveAudV)}')

plt.sca(axFtTeA)
circle8 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveTeA = np.sum(ftSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nSpeechResponsiveTeA - nFTselectiveTeA,  nFTselectiveTeA], colors = [colorNotSelective, colorFtSelective])
axFtTeA.add_artist(circle8)
plt.title(f'TeA,\n n = {int(nSpeechResponsiveTeA)}')
plt.legend(labels = ['FT non-selective', 'FT selective'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))


labelPosX = [0.05, 0.43, 0.63] # Horiz position for panel labels
labelPosY = [0.94, 0.72,  0.46, 0.24]    # Vert position for panel labels

axColorMapVOT.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axVotAudD.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axVotDV.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axVotAP.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axColorMapFT.annotate('E', xy=(labelPosX[0],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axFtAudD.annotate('F', xy=(labelPosX[1],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axFtDV.annotate('G', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axFtAP.annotate('H', xy=(labelPosX[2],labelPosY[3]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
