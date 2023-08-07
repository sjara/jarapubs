"""
This creates figure 2 for 2022paspeech:
 A. Cartoon of headfixed, awake mouse ephys
 B. Diagram of sound matrix
 C. Histology image of recording track
 D. Scatter plot of recording location of each cell. Add AC areas image in inkscape.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import colorpalette as cp
import figparams
import studyparams
from scipy import stats
from importlib import reload

reload(figparams)

FIGNAME = 'selectivityIndices'
SAVE_FIGURE = 1
STATSUMMARY = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'figure_neuropix_methods' # Do not include extension
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


labelPosX = [0.04, 0.28, 0.67] # Horiz position for panel labels
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


if 1:
    # Plot a grid with cell borders and labels
    plt.sca(axSoundMatrix)
    corner_labels = [ ['/da/', '/ta/'], ['/ba/', '/pa/'] ]
    VOTstep = 20
    VOTvals = np.arange(0, 60+VOTstep, VOTstep)
    FTstep = 6
    FTvals = np.arange(-9, 9+FTstep, FTstep)
    for indFT, thisFT in enumerate(FTvals):
        for indVOT, thisVOT in enumerate(VOTvals):
            if indFT in [0, 3] or indVOT in [0, 3]:
                facecolor = '#fce94f'  # Butter2
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
quadrantBoundsDV = figData['quadrantBoundsDV']
quadrantBoundsAP = figData['quadrantBoundsAP']
quadrantBoundsDV_AAtransform = figData['quadrantBoundsDV_AAtransform']
quadrantBoundsAP_AAtransform = figData['quadrantBoundsAP_AAtransform']
quadrantLabels = figData['quadrantLabels']
quadrantTotals = figData['quadrantTotals']
quadrantsSpeechResponsive = figData['quadrantsSpeechResponsive']
quadrantsSoundResponsive = figData['quadrantsSoundResponsive']

APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)


plt.sca(axCellLocs)
nonResp = plt.scatter(z_coords_jittered[~soundResponsive & isCortical], y_coords[~soundResponsive & isCortical], c = colorNotAud, s=6)
soundResp = plt.scatter(z_coords_jittered[soundResponsive & ~speechResponsive & isCortical], y_coords[soundResponsive & ~speechResponsive & isCortical], c = colorSoundResp, s = 6)
speechResp = plt.scatter(z_coords_jittered[speechResponsive & isCortical] , y_coords[speechResponsive & isCortical], c = colorSpeechResp, s=6)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
#plt.legend([nonResp, soundResp, speechResp], ['not sound responsive', 'sound responsive','speech responsive'], loc = 'upper left', markerscale = 2 , bbox_to_anchor = (0, 1.15))
plt.legend([nonResp, soundResp, speechResp], ['Not sound responsive', 'Sound responsive','Speech responsive'], loc = 'upper center', markerscale=2, handletextpad=0.25, bbox_to_anchor = (0.5, 1.07))
axCellLocs.spines["right"].set_visible(False)
axCellLocs.spines["top"].set_visible(False)
axCellLocs.set_aspect('equal')

##-- Donut Plots Atlas-defined Areas--
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0]))
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1]))
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2]) )
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3]))

nSoundResponsiveAudP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[0]))
nSoundResponsiveAudD = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[1]))
nSoundResponsiveAudV = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[2]))
nSoundResponsiveTeA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[3]))

nCellsAudP = np.sum((recordingAreaName == audCtxAreas[0]))
nCellsAudD = np.sum((recordingAreaName == audCtxAreas[1]))
nCellsAudV = np.sum((recordingAreaName == audCtxAreas[2]))
nCellsTeA = np.sum((recordingAreaName == audCtxAreas[3]))


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
nSpeechResponsiveDP = np.sum(quadrantsSpeechResponsive[0])
nSpeechResponsiveDA = np.sum(quadrantsSpeechResponsive[1])
nSpeechResponsiveVP = np.sum(quadrantsSpeechResponsive[2])
nSpeechResponsiveVA = np.sum(quadrantsSpeechResponsive[3])

nSoundResponsiveDP = np.sum(quadrantsSoundResponsive[0])
nSoundResponsiveDA = np.sum(quadrantsSoundResponsive[1])
nSoundResponsiveVP = np.sum(quadrantsSoundResponsive[2])
nSoundResponsiveVA = np.sum(quadrantsSoundResponsive[3])

nCellsDP = quadrantTotals[0]
nCellsDA = quadrantTotals[1]
nCellsVP = quadrantTotals[2]
nCellsVA = quadrantTotals[3]

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
    oddsratio, pvalFracResponsive_AudPvsTea = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveTeA],[nSoundResponsiveAudP - nSpeechResponsiveAudP, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudDvsTea = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveTeA],[nSoundResponsiveAudD - nSpeechResponsiveAudD, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudVvsTea = stats.fisher_exact(np.array([[nSpeechResponsiveAudV, nSpeechResponsiveTeA],[nSoundResponsiveAudV - nSpeechResponsiveAudV, nSoundResponsiveTeA - nSpeechResponsiveTeA]]))

    # Test proportion of speech responsive cells between areas (of all cells)
    oddsratio, pvalFracResponsive_AudPvsAudD_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudD],[nCellsAudP - nSpeechResponsiveAudP, nCellsAudD - nSpeechResponsiveAudD]]))
    oddsratio, pvalFracResponsive_AudPvsAudV_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveAudV],[nCellsAudP - nSpeechResponsiveAudP, nCellsAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudDvsAudV_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveAudV],[nCellsAudD - nSpeechResponsiveAudD, nCellsAudV - nSpeechResponsiveAudV]]))
    oddsratio, pvalFracResponsive_AudPvsTea_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudP, nSpeechResponsiveTeA],[nCellsAudP - nSpeechResponsiveAudP, nCellsTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudDvsTea_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudD, nSpeechResponsiveTeA],[nCellsAudD - nSpeechResponsiveAudD, nCellsTeA - nSpeechResponsiveTeA]]))
    oddsratio, pvalFracResponsive_AudVvsTea_allcells = stats.fisher_exact(np.array([[nSpeechResponsiveAudV, nSpeechResponsiveTeA],[nCellsAudV - nSpeechResponsiveAudV, nCellsTeA - nSpeechResponsiveTeA]]))

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
    print(f'AudP vs Tea p = {np.round(pvalFracResponsive_AudPvsTea,3)}')
    print(f'AudD vs Tea p = {np.round(pvalFracResponsive_AudDvsTea,3)}')
    print(f'AudV vs Tea p = {np.round(pvalFracResponsive_AudVvsTea,3)}')
    print('--Frac speech responsive of total cells --')
    print(f'AudP vs AudD p = {np.round(pvalFracResponsive_AudPvsAudD_allcells,3)}')
    print(f'AudP vs AudV p = {np.round(pvalFracResponsive_AudPvsAudV_allcells,3)}')
    print(f'AudD vs AudV p = {np.round(pvalFracResponsive_AudDvsAudV_allcells,3)}')
    print(f'AudP vs Tea p = {np.round(pvalFracResponsive_AudPvsTea_allcells,3)}')
    print(f'AudD vs Tea p = {np.round(pvalFracResponsive_AudDvsTea_allcells,3)}')
    print(f'AudV vs Tea p = {np.round(pvalFracResponsive_AudVvsTea_allcells,3)}')
    print('--quadrant-defined areas')
    print(f'DP Quadrant: total n = {len(quadrantsSoundResponsive[0])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[0])} ({np.round((np.sum(quadrantsSoundResponsive[0])/len(quadrantsSoundResponsive[0]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[0])} ({np.round((np.sum(quadrantsSpeechResponsive[0])/np.sum(quadrantsSoundResponsive[0]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[0])/len(quadrantsSoundResponsive[0]))*100,1)}% of total)')
    print(f'DA Quadrant: total n = {len(quadrantsSoundResponsive[1])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[1])} ({np.round((np.sum(quadrantsSoundResponsive[1])/len(quadrantsSoundResponsive[1]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[1])} ({np.round((np.sum(quadrantsSpeechResponsive[1])/np.sum(quadrantsSoundResponsive[1]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[1])/len(quadrantsSoundResponsive[1]))*100,1)}% of total)')
    print(f'VP Quadrant: total n = {len(quadrantsSoundResponsive[2])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[2])} ({np.round((np.sum(quadrantsSoundResponsive[2])/len(quadrantsSoundResponsive[2]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[2])} ({np.round((np.sum(quadrantsSpeechResponsive[2])/np.sum(quadrantsSoundResponsive[2]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[2])/len(quadrantsSoundResponsive[2]))*100,1)}% of total)')
    print(f'VA Quadrant: total n = {len(quadrantsSoundResponsive[3])}, n soundResponsive = {np.sum(quadrantsSoundResponsive[3])} ({np.round((np.sum(quadrantsSoundResponsive[3])/len(quadrantsSoundResponsive[3]))*100,1)}%), n speech responsive = {np.sum(quadrantsSpeechResponsive[3])} ({np.round((np.sum(quadrantsSpeechResponsive[3])/np.sum(quadrantsSoundResponsive[3]))*100,1)}% of sound responsive, {np.round((np.sum(quadrantsSpeechResponsive[3])/len(quadrantsSoundResponsive[3]))*100,1)}% of total)')
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
