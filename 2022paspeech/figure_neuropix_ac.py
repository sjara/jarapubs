'''
This creates figure 2 for 2022paspeech:
 A. Cartoon of headfixed, awake mouse ephys
 B. Diagram of sound matrix
 C. Histology image of recording track
 D. Scatter plot of recording location of each cell. Add AC areas image in inkscape.
'''
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
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.font_manager
matplotlib.font_manager.fontManager.addfont('C:\\Users\\jenny\\anaconda3\\Lib\\site-packages\\matplotlib\\mpl-data\\fonts\\ttf\\Helvetica.ttf')
reload(figparams)

FIGNAME = 'selectivityIndices'
SAVE_FIGURE = 0
STATSUMMARY = 0
outputDir = 'C:/Users/jenny/tmp/'
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


labelPosX = [0.07, 0.28, 0.58] # Horiz position for panel labels
labelPosY = [0.93, 0.75, 0.375]    # Vert position for panel labels

gsMain = gridspec.GridSpec(2, 1, height_ratios=[0.15, 0.85])
gsMain.update(left=0.1, right=0.97, top=0.9, bottom=0.08, wspace=0.3, hspace=0.3)
gsTop = gsMain[0].subgridspec(1,3, width_ratios = [0.2, 0.35, 0.45])

axCartoon = plt.subplot(gsTop[0,0])
axCartoon.set_axis_off()

axHist = plt.subplot(gsTop[0,1])
axHist.set_axis_off()

axSoundMatrix = plt.subplot(gsTop[0,2])
axSoundMatrix.set_axis_off()

gsBottom = gsMain[1].subgridspec(4,3, width_ratios = [0.6, 0.2, 0.2])

axCellLocs = plt.subplot(gsBottom[:,0])

axDonutAudP = plt.subplot(gsBottom[0,1])
axDonutAudD = plt.subplot(gsBottom[0,2])
axDonutAudV = plt.subplot(gsBottom[1,1])
axDonutTeA = plt.subplot(gsBottom[1,2])


axDonutDP = plt.subplot(gsBottom[2,1])
axDonutDA = plt.subplot(gsBottom[2,2])
axDonutVP = plt.subplot(gsBottom[3,1])
axDonutVA = plt.subplot(gsBottom[3,2])


axCartoon.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axSoundMatrix.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axHist.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axCellLocs.annotate('D', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axDonutAudP.annotate('E', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axDonutDP.annotate('F', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

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
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.legend([nonResp, soundResp, speechResp], ['not sound responsive', 'sound responsive','speech responsive'], loc = 'upper left', markerscale = 2 , bbox_to_anchor = (0, 1.15))
axCellLocs.spines["right"].set_visible(False)
axCellLocs.spines["top"].set_visible(False)

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


plt.sca(axDonutAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudP, nSoundResponsiveAudP - nSpeechResponsiveAudP, nCellsAudP - nSoundResponsiveAudP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudP.add_artist(circle1)
plt.title(f'AudP,\n n = {nCellsAudP}', pad = -0.3 , fontsize = fontSizeTicks)


plt.sca(axDonutAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudD, nSoundResponsiveAudD - nSpeechResponsiveAudD, nCellsAudD - nSoundResponsiveAudD], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudD.add_artist(circle2)
plt.title(f'AudD,\n n = {nCellsAudD}', pad = -0.3, fontsize = fontSizeTicks)


plt.sca(axDonutAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveAudV, nSoundResponsiveAudV - nSpeechResponsiveAudV, nCellsAudV - nSoundResponsiveAudV], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutAudV.add_artist(circle3)
plt.title(f'AudV,\n n = {nCellsAudV}', pad = -0.3, fontsize = fontSizeTicks)

plt.sca(axDonutTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveTeA, nSoundResponsiveTeA - nSpeechResponsiveTeA, nCellsTeA - nSoundResponsiveTeA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutTeA.add_artist(circle4)
plt.title(f'TeA,\n n = {nCellsTeA}', pad = -0.3, fontsize = fontSizeTicks)


##-- Donut Plots Quartiles--
nSpeechResponsiveDP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0]))
nSpeechResponsiveDA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1]))
nSpeechResponsiveVP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2]) )
nSpeechResponsiveVA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3]))

nSoundResponsiveDP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[0]))
nSoundResponsiveDA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[1]))
nSoundResponsiveVP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[2]))
nSoundResponsiveVA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[3]))

nCellsDP = np.sum((recordingAreaName == audCtxAreas[0]))
nCellsDA = np.sum((recordingAreaName == audCtxAreas[1]))
nCellsVP = np.sum((recordingAreaName == audCtxAreas[2]))
nCellsVA = np.sum((recordingAreaName == audCtxAreas[3]))


plt.sca(axDonutDP)
circle5 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveDP, nSoundResponsiveDP - nSpeechResponsiveDP, nCellsDP - nSoundResponsiveDP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutDP.add_artist(circle5)
plt.title(f'DP,\n n = {nCellsDP}', pad = -0.3 , fontsize = fontSizeTicks)


plt.sca(axDonutDA)
circle6 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveDA, nSoundResponsiveDA - nSpeechResponsiveDA, nCellsDA - nSoundResponsiveDA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutDA.add_artist(circle6)
plt.title(f'DA,\n n = {nCellsDA}', pad = -0.3, fontsize = fontSizeTicks)


plt.sca(axDonutVP)
circle7 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveVP, nSoundResponsiveVP - nSpeechResponsiveVP, nCellsVP - nSoundResponsiveVP], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutVP.add_artist(circle7)
plt.title(f'VP,\n n = {nCellsVP}', pad = -0.3, fontsize = fontSizeTicks)

plt.sca(axDonutVA)
circle8 = plt.Circle((0,0), 0.7, color = 'white')
plt.pie([nSpeechResponsiveVA, nSoundResponsiveVA - nSpeechResponsiveVA, nCellsVA - nSoundResponsiveVA], colors = [colorSpeechResp, colorSoundResp, colorNotAud])
axDonutVA.add_artist(circle8)
plt.title(f'VA,\n n = {nCellsVA}', pad = -0.3, fontsize = fontSizeTicks)



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


    print('--Stats Summary--')
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
