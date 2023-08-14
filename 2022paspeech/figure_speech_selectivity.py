"""
This creates figure 1 for 2022paspeech:
 A. plot of VOT selectivity indices in A-P/D-V space. Add AC areas image in inkscape
 B. Donut plots of fraction of VOT selective cells among speech responsive cells for each Aud Area
 C. VOT selectivity indices binned in D-V space
 D. VOT selectivity indices binned in A-P space
 E-H: as in A-D, for FT
"""

import sys
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


FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
shuffledDataFile = 'data_shuffledSIs.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
SAVE_FIGURE = 1
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'plots_speech_selectivity' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7.5, 9] # In inches
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
subjects = np.unique(figData['subject'])
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
isCortical = figData['isCortical']
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
quadrantBoundsDV_AApix = figData['quadrantBoundsDV_AApix']
quadrantBoundsAP_AApix = figData['quadrantBoundsAP_AApix']
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



plt.clf()
gsMain = gridspec.GridSpec(2, 1)
gsMain.update(left=0.08, right=0.98, top=0.99, bottom=0.04, wspace=0.25, hspace=0.2)
#plt.subplots_adjust(top = 0.9, bottom = 0.05, hspace = 0.45, left = 0.05)

gsVOT = gsMain[0].subgridspec(4, 3, width_ratios=[0.55, 0.15, 0.3], wspace=0.55)
axColorMapVOT = plt.subplot(gsVOT[:,0])
axDVAPmapVot = gsVOT[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45], hspace=0.4)
axQuadsVot = plt.subplot(axDVAPmapVot[0,0])
axByAnimalVot = plt.subplot(axDVAPmapVot[1,0])
axVotDonuts = gsVOT[0:,1].subgridspec(4, 1, hspace=0.2)
axVotAudP = plt.subplot(axVotDonuts[1,0])
axVotAudD = plt.subplot(axVotDonuts[0,0])
axVotAudV = plt.subplot(axVotDonuts[2,0])
axVotTeA = plt.subplot(axVotDonuts[3,0])
# -- Move Donuts a little to the left and closer together --
xoffset = -0.02
yoffset = 0.01
for indax, oneax in enumerate([axVotAudD, axVotAudP, axVotAudV, axVotTeA]):
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0+xoffset, axpos.y0+(indax-1)*yoffset, axpos.width, axpos.height])

#gsFT = gsMain[1].subgridspec(4, 3, width_ratios=[0.5, 0.15, 0.35], wspace=0.55)
gsFT = gsMain[1].subgridspec(4, 3, width_ratios=[0.55, 0.15, 0.3], wspace=0.55)
axColorMapFT = plt.subplot(gsFT[:,0])
axDVAPmapFt = gsFT[:,2].subgridspec(2,1, height_ratios = [0.55, 0.45], hspace=0.4)
axQuadsFt = plt.subplot(axDVAPmapFt[0,0])
axByAnimalFt = plt.subplot(axDVAPmapFt[1,0])
axFtDonuts = gsFT[0:,1].subgridspec(4, 1, hspace=0.2)
axFtAudP = plt.subplot(axFtDonuts[1,0])
axFtAudD = plt.subplot(axFtDonuts[0,0])
axFtAudV = plt.subplot(axFtDonuts[2,0])
axFtTeA = plt.subplot(axFtDonuts[3,0])
# -- Move Donuts a little to the left and closer together --
xoffset = -0.02
yoffset = 0.01
for indax, oneax in enumerate([axFtAudD, axFtAudP, axFtAudV, axFtTeA]):
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0+xoffset, axpos.y0+(indax-1)*yoffset, axpos.width, axpos.height])

# -- Move third column a little higher --
yoffset = 0.01
for indax, oneax in enumerate([axQuadsVot, axByAnimalVot, axQuadsFt, axByAnimalFt]):
    axpos = oneax.get_position()
    oneax.set_position([axpos.x0, axpos.y0+yoffset, axpos.width, axpos.height])


votSelective = (pvalPermutationtestVot < 0.05) & speechResponsive
ftSelective = (pvalPermutationtestFt < 0.05) & speechResponsive
singleSelective = np.logical_xor(votSelective, ftSelective)
mixedSelective = votSelective & ftSelective


#plt.show();    sys.exit()



#APtickLocs = np.array([176, 196, 216])
APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
#DVtickLocs = np.array([190, 140, 90])
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
#plt.suptitle('VOT selectivity by location', fontsize = fontSizeTitles)

plt.sca(axColorMapVOT)
plt.scatter(z_coords_jittered[speechResponsive & isCortical], y_coords[speechResponsive & isCortical],
            c=bestSelectivityIndexVot[speechResponsive & isCortical],
            cmap=newMapVOT, s=3, vmin=0, vmax=1)
#plt.scatter(z_coords_jittered, y_coords, c = bestSelectivityIndexVot, cmap = newMapVOT, s = 3)
#plt.ylim(215,60)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
#plt.xlim(165,225)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)
#cbar = plt.colorbar(shrink=0.8, ticks=[0.2, 0.4, 0.6, 0.8, 1.0])
cbar = plt.colorbar(ticks=[0.2, 0.4, 0.6, 0.8, 1.0], shrink=0.8, pad=0.1)
cbar.set_label('VOT Selectivity Index', rotation=270, labelpad=12)
htitle = plt.title('VOT selectivity', fontsize = fontSizeTitles, y=0.95)
htitle.set_position(htitle.get_position()+np.array([0.08, 0]))
axColorMapVOT.spines["right"].set_visible(False)
axColorMapVOT.spines["top"].set_visible(False)
axColorMapVOT.set_aspect('equal')

titleY = 0.9

plt.sca(axVotAudP)
circle1 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudP = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[0] ))
nSoundResponsiveAudP = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[0]))
nTotalAudP = np.sum(recordingAreaName == audCtxAreas[0])
nVOTselectiveAudP = np.sum(votSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nVOTselectiveAudP, nTotalAudP - nVOTselectiveAudP], colors = [colorVotSelective, colorNotSelective])
axVotAudP.add_artist(circle1)
#plt.title(f'AudP,\n n = {int(nTotalAudP)}')
plt.title(f'AudP (n={int(nTotalAudP)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axVotAudD)
circle2 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudD = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[1] ))
nSoundResponsiveAudD = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[1]))
nTotalAudD = np.sum(recordingAreaName == audCtxAreas[1])
nVOTselectiveAudD = np.sum(votSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([ nVOTselectiveAudD, nTotalAudD - nVOTselectiveAudD], colors = [colorVotSelective, colorNotSelective])
axVotAudD.add_artist(circle2)
#plt.title(f'AudD,\n n = {int(nTotalAudD)}')
plt.title(f'AudD (n={int(nTotalAudD)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axVotAudV)
circle3 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveAudV = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[2] ))
nSoundResponsiveAudV = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[2]))
nTotalAudV = np.sum(recordingAreaName == audCtxAreas[2])
nVOTselectiveAudV = np.sum(votSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nVOTselectiveAudV, nTotalAudV - nVOTselectiveAudV], colors = [colorVotSelective, colorNotSelective])
axVotAudV.add_artist(circle3)
#plt.title(f'AudV,\n n = {int(nTotalAudV)}')
plt.title(f'AudV (n={int(nTotalAudV)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axVotTeA)
circle4 = plt.Circle((0,0), 0.7, color = 'white')
nSpeechResponsiveTeA = np.sum(speechResponsive & (recordingAreaName == audCtxAreas[3] ))
nSoundResponsiveTeA = np.sum(soundResponsive & (recordingAreaName == audCtxAreas[3] ))
nTotalTeA = np.sum(recordingAreaName == audCtxAreas[3])
nVOTselectiveTeA = np.sum(votSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nVOTselectiveTeA, nTotalTeA - nVOTselectiveTeA], colors = [colorVotSelective, colorNotSelective])
axVotTeA.add_artist(circle4)
#plt.title(f'TeA,\n n = {int(nTotalTeA)}')
plt.title(f'TeA (n={int(nTotalTeA)})', y=titleY, fontsize=fontSizeLabels)
#plt.legend(labels = ['VOT selective', 'VOT non-sel.'], loc = 'lower center', bbox_to_anchor = (0.1, -0.9))
plt.legend(labels = ['VOT selective', 'VOT non-sel.'], loc='lower center', handlelength=1.5,
           bbox_to_anchor = (0.5, -0.75))

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
plt.annotate(f'DP\n{np.round(fracVOTselectiveDP*100,1)}%\nn = {nTotalDP}', (0.25, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'w', fontweight = 'normal')
plt.annotate(f'DA\n{np.round(fracVOTselectiveDA*100,1)}%\nn = {nTotalDA}', (0.75, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'w', fontweight = 'normal')
plt.annotate(f'VP\n{np.round(fracVOTselectiveVP*100,1)}%\nn = {nTotalVP}', (0.25, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k', fontweight = 'normal')
plt.annotate(f'VA\n{np.round(fracVOTselectiveVA*100,1)}%\nn = {nTotalVA}', (0.75, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k', fontweight = 'normal')
axQuadsVot.spines["right"].set_visible(False)
axQuadsVot.spines["top"].set_visible(False)
axQuadsVot.set_aspect('equal')
plt.xticks([0,0.5,1], labels = np.round(quadrantBoundsAP, 2))
plt.yticks([0,0.5,1], labels = np.round(quadrantBoundsDV[::-1], 2))
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)


plt.sca(axByAnimalVot)
fracVotSelectiveByAnimal = quadrantsVotSelectiveByAnimal/quadrantTotalsByAnimal
fracVotSelectiveByAnimal[quadrantTotalsByAnimal<=2] = np.nan
medianVotSelectiveByAnimal = np.nanmedian(fracVotSelectiveByAnimal, axis = 0)
quadXcoords = np.array([1,2,3,4])
plt.bar(quadXcoords, [medianVotSelectiveByAnimal[0], medianVotSelectiveByAnimal[1], medianVotSelectiveByAnimal[2], medianVotSelectiveByAnimal[3]], facecolor = colorVotSelective)
plt.xticks(quadXcoords, ['DP', 'DA', 'VP', 'VA'])
plt.ylabel('Fraction VOT selective')
plt.xlabel('AC regions')
for indAnimal, thisAnimal in enumerate(subjects):
    plt.plot(quadXcoords[quadrantTotalsByAnimal[indAnimal]>2], fracVotSelectiveByAnimal[indAnimal,quadrantTotalsByAnimal[indAnimal]>2], c = colorNotSelective, marker = 'o', ms = 4, alpha = 0.5)
axByAnimalVot.spines["right"].set_visible(False)
axByAnimalVot.spines["top"].set_visible(False)


plt.sca(axColorMapFT)
#plt.scatter(z_coords_jittered[speechResponsive & isCortical], y_coords[speechResponsive& isCortical], c = bestSelectivityIndexFt[speechResponsive& isCortical], cmap = newMapFT, s = 3)
plt.scatter(z_coords_jittered[speechResponsive & isCortical], y_coords[speechResponsive& isCortical],
            c=bestSelectivityIndexFt[speechResponsive& isCortical],
            cmap=newMapFT, s=3, vmin=0, vmax=1)
#plt.scatter(z_coords_jittered, y_coords, c = bestSelectivityIndexFt, cmap = newMapFT, s = 3)
#plt.ylim(215,60)
#plt.xlim(165,225)
plt.ylim(220,40)
plt.xlim(146, 246)
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)
plt.xticks(APtickLocs, APtickLabels)
plt.yticks(DVtickLocs, DVtickLabels)
#cbar = plt.colorbar(shrink=0.8, ticks=[0.2, 0.4, 0.6, 0.8, 1.0], extend='max', extendrect=True, extendfrac=0.22)
cbar = plt.colorbar(ticks=[0.2, 0.4, 0.6, 0.8, 1.0], shrink=0.8, pad=0.1)
cbar.set_label('FT Selectivity Index', rotation=270, labelpad=12)
htitle = plt.title('FT selectivity', fontsize = fontSizeTitles, y=0.95)
htitle.set_position(htitle.get_position()+np.array([0.08, 0]))
axColorMapFT.set_aspect('equal')
axColorMapFT.spines["right"].set_visible(False)
axColorMapFT.spines["top"].set_visible(False)



plt.sca(axFtAudP)
circle5 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudP = np.sum(ftSelective[recordingAreaName == audCtxAreas[0]])
plt.pie([nFTselectiveAudP, nTotalAudP - nFTselectiveAudP], colors = [colorFtSelective, colorNotSelective])
axFtAudP.add_artist(circle5)
plt.title(f'AudP (n={int(nTotalAudP)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axFtAudD)
circle6 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudD = np.sum(ftSelective[recordingAreaName == audCtxAreas[1]])
plt.pie([nFTselectiveAudD, nTotalAudD - nFTselectiveAudD], colors = [colorFtSelective, colorNotSelective])
axFtAudD.add_artist(circle6)
plt.title(f'AudD (n={int(nTotalAudD)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axFtAudV)
circle7 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveAudV = np.sum(ftSelective[recordingAreaName == audCtxAreas[2]])
plt.pie([nFTselectiveAudV, nTotalAudV - nFTselectiveAudV], colors = [colorFtSelective, colorNotSelective])
axFtAudV.add_artist(circle7)
plt.title(f'AudV (n={int(nTotalAudV)})', y=titleY, fontsize=fontSizeLabels)

plt.sca(axFtTeA)
circle8 = plt.Circle((0,0), 0.7, color = 'white')
nFTselectiveTeA = np.sum(ftSelective[recordingAreaName == audCtxAreas[3]])
plt.pie([nFTselectiveTeA, nTotalTeA - nFTselectiveTeA], colors = [colorFtSelective, colorNotSelective])
axFtTeA.add_artist(circle8)
plt.title(f'TeA (n={int(nTotalTeA)})', y=titleY, fontsize=fontSizeLabels)
plt.legend(labels=['FT selective', 'FT non-sel.'], loc='lower center', handlelength=1.5,
           bbox_to_anchor=(0.5, -0.75))


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
plt.annotate(f'DP\n{np.round(fracFTselectiveDP*100,1)}%\nn = {nTotalDP}', (0.25, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'w')
plt.annotate(f'DA\n{np.round(fracFTselectiveDA*100,1)}%\nn = {nTotalDA}', (0.75, 0.6), fontsize = fontSizeLabels, ha = 'center', color = 'k')
plt.annotate(f'VP\n{np.round(fracFTselectiveVP*100,1)}%\nn = {nTotalVP}', (0.25, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k')
plt.annotate(f'VA\n{np.round(fracFTselectiveVA*100,1)}%\nn = {nTotalVA}', (0.75, 0.1), fontsize = fontSizeLabels, ha = 'center', color = 'k')
axQuadsFt.spines["right"].set_visible(False)
axQuadsFt.spines["top"].set_visible(False)
axQuadsFt.set_aspect('equal')
plt.xticks([0,0.5,1], labels = np.round(quadrantBoundsAP, 2))
plt.yticks([0,0.5,1], labels = np.round(quadrantBoundsDV[::-1], 2))
plt.xlim([-0.1,1.1])
plt.ylim([-0.1,1.1])
plt.ylabel('Ventral-Dorsal (mm)', fontsize = fontSizeLabels)
plt.xlabel('Posterior-Anterior (mm)', fontsize = fontSizeLabels)
#plt.axis('off')


plt.sca(axByAnimalFt)
fracFtSelectiveByAnimal = quadrantsFtSelectiveByAnimal/quadrantTotalsByAnimal
fracFtSelectiveByAnimal[quadrantTotalsByAnimal<=2] = np.nan
medianFTselectiveByAnimal = np.nanmedian(fracFtSelectiveByAnimal, axis = 0)
plt.bar(quadXcoords, [medianFTselectiveByAnimal[0], medianFTselectiveByAnimal[1], medianFTselectiveByAnimal[2], medianFTselectiveByAnimal[3]], facecolor = colorFtSelective)
#plt.bar([1, 2, 3, 4], [fracFTselectiveDP, fracFTselectiveDA, fracFTselectiveVA, fracFTselectiveVP], facecolor = colorFtSelective)
plt.xticks(quadXcoords, ['DP', 'DA', 'VP', 'VA'])
plt.ylabel('Fraction FT selective', fontsize = fontSizeLabels)
plt.xlabel('AC regions', fontsize = fontSizeLabels)
#plt.plot([1,2,3,4], [fracFtSelectiveByAnimal[:,0], fracFtSelectiveByAnimal[:,1], fracFtSelectiveByAnimal[:,3], fracFtSelectiveByAnimal[:,2]], c = colorNotSelective, alpha = 0.5)

for indAnimal, thisAnimal in enumerate(subjects):
    plt.plot(quadXcoords[quadrantTotalsByAnimal[indAnimal]>2], fracFtSelectiveByAnimal[indAnimal,quadrantTotalsByAnimal[indAnimal]>2], c = colorNotSelective, marker = 'o', ms = 4, alpha = 0.5)
axByAnimalFt.spines["right"].set_visible(False)
axByAnimalFt.spines["top"].set_visible(False)



labelPosX = [0.02, 0.50, 0.70] # Horiz position for panel labels
labelPosY = [0.98, 0.74,  0.47, 0.22]    # Vert position for panel labels

axColorMapVOT.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction',
                       fontsize=fontSizePanel, fontweight='bold')
axVotAudD.annotate('B', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axQuadsVot.annotate('C', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction',
                    fontsize=fontSizePanel, fontweight='bold')
axByAnimalVot.annotate('D', xy=(labelPosX[2],labelPosY[1]), xycoords='figure fraction',
                       fontsize=fontSizePanel, fontweight='bold')
axColorMapFT.annotate('E', xy=(labelPosX[0],labelPosY[2]), xycoords='figure fraction',
                      fontsize=fontSizePanel, fontweight='bold')
axFtAudD.annotate('F', xy=(labelPosX[1],labelPosY[2]), xycoords='figure fraction',
                  fontsize=fontSizePanel, fontweight='bold')
axQuadsFt.annotate('G', xy=(labelPosX[2],labelPosY[2]), xycoords='figure fraction',
                   fontsize=fontSizePanel, fontweight='bold')
axByAnimalFt.annotate('H', xy=(labelPosX[2],labelPosY[3]), xycoords='figure fraction',
                      fontsize=fontSizePanel, fontweight='bold')

plt.show()


if STATSUMMARY:
    # -- ATLAS AREAS
    ## VOT selective out of all cells
    oddsratio, pvalFracVotSelective_AudPvsAudD = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveAudD],[nTotalAudP - nVOTselectiveAudP, nTotalAudD - nVOTselectiveAudD]]))
    oddsratio, pvalFracVotSelective_AudPvsAudV = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveAudV],[nTotalAudP - nVOTselectiveAudP, nTotalAudV - nVOTselectiveAudV]]))
    oddsratio, pvalFracVotSelective_AudDvsAudV = stats.fisher_exact(np.array([[nVOTselectiveAudD, nVOTselectiveAudV],[nTotalAudD - nVOTselectiveAudD, nTotalAudV - nVOTselectiveAudV]]))
    oddsratio, pvalFracVotSelective_AudPvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveTeA],[nTotalAudP - nVOTselectiveAudP, nTotalTeA - nVOTselectiveTeA]]))
    oddsratio, pvalFracVotSelective_AudDvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudD, nVOTselectiveTeA],[nTotalAudD - nVOTselectiveAudD, nTotalTeA - nVOTselectiveTeA]]))
    oddsratio, pvalFracVotSelective_AudVvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudV, nVOTselectiveTeA],[nTotalAudV - nVOTselectiveAudV, nTotalTeA - nVOTselectiveTeA]]))

    ## VOT selective out of sound responsive
    oddsratio, pvalFracVotSelective_soundResponsiveAudPvsAudD = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveAudD],[nSoundResponsiveAudP - nVOTselectiveAudP, nSoundResponsiveAudD - nVOTselectiveAudD]]))
    oddsratio, pvalFracVotSelective_soundResponsiveAudPvsAudV = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveAudV],[nSoundResponsiveAudP - nVOTselectiveAudP, nSoundResponsiveAudV - nVOTselectiveAudV]]))
    oddsratio, pvalFracVotSelective_soundResponsiveAudDvsAudV = stats.fisher_exact(np.array([[nVOTselectiveAudD, nVOTselectiveAudV],[nSoundResponsiveAudD - nVOTselectiveAudD, nSoundResponsiveAudV - nVOTselectiveAudV]]))
    oddsratio, pvalFracVotSelective_soundResponsiveAudPvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudP, nVOTselectiveTeA],[nSoundResponsiveAudP - nVOTselectiveAudP, nSoundResponsiveTeA - nVOTselectiveTeA]]))
    oddsratio, pvalFracVotSelective_soundResponsiveAudDvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudD, nVOTselectiveTeA],[nSoundResponsiveAudD - nVOTselectiveAudD, nSoundResponsiveTeA - nVOTselectiveTeA]]))
    oddsratio, pvalFracVotSelective_soundResponsiveAudVvsTeA = stats.fisher_exact(np.array([[nVOTselectiveAudV, nVOTselectiveTeA],[nSoundResponsiveAudV - nVOTselectiveAudV, nSoundResponsiveTeA - nVOTselectiveTeA]]))

    ## FT selective out of all cells
    oddsratio, pvalFracFtSelective_AudPvsAudD = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveAudD],[nTotalAudP - nFTselectiveAudP, nTotalAudD - nFTselectiveAudD]]))
    oddsratio, pvalFracFtSelective_AudPvsAudV = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveAudV],[nTotalAudP - nFTselectiveAudP, nTotalAudV - nFTselectiveAudV]]))
    oddsratio, pvalFracFtSelective_AudDvsAudV = stats.fisher_exact(np.array([[nFTselectiveAudD, nFTselectiveAudV],[nTotalAudD - nFTselectiveAudD, nTotalAudV - nFTselectiveAudV]]))
    oddsratio, pvalFracFtSelective_AudPvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveTeA],[nTotalAudP - nFTselectiveAudP, nTotalTeA - nFTselectiveTeA]]))
    oddsratio, pvalFracFtSelective_AudDvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudD, nFTselectiveTeA],[nTotalAudD - nFTselectiveAudD, nTotalTeA - nFTselectiveTeA]]))
    oddsratio, pvalFracFtSelective_AudVvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudV, nFTselectiveTeA],[nTotalAudV - nFTselectiveAudV, nTotalTeA - nFTselectiveTeA]]))

    ## FT selective out of sound responsive
    oddsratio, pvalFracFtSelective_soundResponsiveAudPvsAudD = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveAudD],[nSoundResponsiveAudP - nFTselectiveAudP, nSoundResponsiveAudD - nFTselectiveAudD]]))
    oddsratio, pvalFracFtSelective_soundResponsiveAudPvsAudV = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveAudV],[nSoundResponsiveAudP - nFTselectiveAudP, nSoundResponsiveAudV - nFTselectiveAudV]]))
    oddsratio, pvalFracFtSelective_soundResponsiveAudDvsAudV = stats.fisher_exact(np.array([[nFTselectiveAudD, nFTselectiveAudV],[nSoundResponsiveAudD - nFTselectiveAudD, nSoundResponsiveAudV - nFTselectiveAudV]]))
    oddsratio, pvalFracFtSelective_soundResponsiveAudPvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudP, nFTselectiveTeA],[nSoundResponsiveAudP - nFTselectiveAudP, nSoundResponsiveTeA - nFTselectiveTeA]]))
    oddsratio, pvalFracFtSelective_soundResponsiveAudDvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudD, nFTselectiveTeA],[nSoundResponsiveAudD - nFTselectiveAudD, nSoundResponsiveTeA - nFTselectiveTeA]]))
    oddsratio, pvalFracFtSelective_soundResponsiveAudVvsTeA = stats.fisher_exact(np.array([[nFTselectiveAudV, nFTselectiveTeA],[nSoundResponsiveAudV - nFTselectiveAudV, nSoundResponsiveTeA - nFTselectiveTeA]]))





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
        quadrantComparFt = np.ones(nQuadCompar)
        a = -1
        for indBin, thisBin in enumerate(quadrantsFT):
            nBinCompar = 4 - indBin -1
            for x in range(nBinCompar):
                a = a+1
                ustat, pvalMannU = stats.mannwhitneyu(quadrantsFT[indBin], quadrantsFT[x + indBin + 1])
                quadrantComparFt[a] = pvalMannU

    ##--Test frac VOT selective of total cells
    quadrantComparFracVotSelective_allcells = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsVotSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsVotSelective[x + indBin + 1])],[np.sum(quadrantTotals[indBin]) -np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantTotals[x + indBin + 1]) - np.sum(quadrantsVotSelective[x + indBin + 1])]]))
            quadrantComparFracVotSelective_allcells[a] = pvalFracVOTSelective

    ##--Test frac FT selective of total cells
    quadrantComparFracFtSelective_allcells = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsFtSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsFtSelective[x + indBin + 1])],[np.sum(quadrantTotals[indBin]) - np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantTotals[x + indBin + 1]) - np.sum(quadrantsFtSelective[x + indBin + 1])]]))
            quadrantComparFracFtSelective_allcells[a] = pvalFracFTSelective

    ##--Test frac VOT selective of sound responsive
    quadrantComparFracVotSelective_soundResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsVotSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsVotSelective[x + indBin + 1])],[np.sum(quadrantsSoundResponsive[indBin]) -np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsVotSelective[x + indBin + 1])]]))
            quadrantComparFracVotSelective_soundResponsive[a] = pvalFracVOTSelective

    ##--Test frac FT selective of sound responsive
    quadrantComparFracFtSelective_soundResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsFtSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsFtSelective[x + indBin + 1])],[np.sum(quadrantsSoundResponsive[indBin]) - np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsSoundResponsive[x + indBin + 1]) - np.sum(quadrantsFtSelective[x + indBin + 1])]]))
            quadrantComparFracFtSelective_soundResponsive[a] = pvalFracFTSelective

    ##--Test frac VOT selective of speech responsive
    quadrantComparFracVotSelective_speechResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsVotSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracVOTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsVotSelective[x + indBin + 1])],[np.sum(quadrantsSpeechResponsive[indBin]) -np.sum(quadrantsVotSelective[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1]) - np.sum(quadrantsVotSelective[x + indBin + 1])]]))
            quadrantComparFracVotSelective_speechResponsive[a] = pvalFracVOTSelective

    ##--Test frac FT selective of speech responsive
    quadrantComparFracFtSelective_speechResponsive = np.ones(nQuadCompar)
    a = -1
    for indBin, thisBin in enumerate(quadrantsFtSelective):
        nBinCompar = 4 - indBin -1
        for x in range(nBinCompar):
            a = a+1
            oddsratio, pvalFracFTSelective = stats.fisher_exact(np.array([[np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsFtSelective[x + indBin + 1])],[np.sum(quadrantsSpeechResponsive[indBin]) - np.sum(quadrantsFtSelective[indBin]), np.sum(quadrantsSpeechResponsive[x + indBin + 1]) - np.sum(quadrantsFtSelective[x + indBin + 1])]]))
            quadrantComparFracFtSelective_speechResponsive[a] = pvalFracFTSelective


        '''
        plt.sca(axFtAudP)
        circle5 = plt.Circle((0,0), 0.7, color = 'white')
        nFTselectiveAudP = np.sum(ftSelective[recordingAreaName == audCtxAreas[0]])
        plt.pie([nFTselectiveAudP, nTotalAudP - nFTselectiveAudP], colors = [colorFtSelective, colorNotSelective])
        axFtAudP.add_artist(circle5)
        plt.title(f'AudP,\n n = {int(nTotalAudP)}')
        '''

    print('--Stat Summary --')
    print('--atlas-defined areas--')
    print(f'AudP: nTotal = {nTotalAudP}, n Sound responsive = {nSoundResponsiveAudP}, n Speech responsive = {nSpeechResponsiveAudP}, n VOT selective = {nVOTselectiveAudP}, ({np.round((nVOTselectiveAudP/nTotalAudP)*100, 1)}% of total cells, {np.round((nVOTselectiveAudP/nSoundResponsiveAudP)*100, 1)}% of sound responsive, {np.round((nVOTselectiveAudP/nSpeechResponsiveAudP)*100, 1)}% of speech responsive), n FT selective = {nFTselectiveAudP}, ({np.round((nFTselectiveAudP/nTotalAudP)*100, 1)}% of total cells, {np.round((nFTselectiveAudP/nSoundResponsiveAudP)*100, 1)}% of sound responsive, {np.round((nFTselectiveAudP/nSpeechResponsiveAudP)*100, 1)}% of speech responsive)')
    print(f'AudD: nTotal = {nTotalAudD}, n Sound responsive = {nSoundResponsiveAudD}, n Speech responsive = {nSpeechResponsiveAudD}, n VOT selective = {nVOTselectiveAudD}, ({np.round((nVOTselectiveAudD/nTotalAudD)*100, 1)}% of total cells, {np.round((nVOTselectiveAudD/nSoundResponsiveAudD)*100, 1)}% of sound responsive, {np.round((nVOTselectiveAudD/nSpeechResponsiveAudD)*100, 1)}% of speech responsive), n FT selective = {nFTselectiveAudD}, ({np.round((nFTselectiveAudD/nTotalAudD)*100, 1)}% of total cells, {np.round((nFTselectiveAudD/nSoundResponsiveAudD)*100, 1)}% of sound responsive, {np.round((nFTselectiveAudD/nSpeechResponsiveAudD)*100, 1)}% of speech responsive)')
    print(f'AudV: nTotal = {nTotalAudV}, n Sound responsive = {nSoundResponsiveAudV}, n Speech responsive = {nSpeechResponsiveAudV}, n VOT selective = {nVOTselectiveAudV}, ({np.round((nVOTselectiveAudV/nTotalAudV)*100, 1)}% of total cells, {np.round((nVOTselectiveAudV/nSoundResponsiveAudV)*100, 1)}% of sound responsive, {np.round((nVOTselectiveAudV/nSpeechResponsiveAudV)*100, 1)}% of speech responsive), n FT selective = {nFTselectiveAudV}, ({np.round((nFTselectiveAudV/nTotalAudV)*100, 1)}% of total cells, {np.round((nFTselectiveAudV/nSoundResponsiveAudV)*100, 1)}% of sound responsive, {np.round((nFTselectiveAudV/nSpeechResponsiveAudV)*100, 1)}% of speech responsive)')
    print(f'TeA: nTotal = {nTotalTeA}, n Sound responsive = {nSoundResponsiveTeA}, n Speech responsive = {nSpeechResponsiveTeA}, n VOT selective = {nVOTselectiveTeA}, ({np.round((nVOTselectiveTeA/nTotalTeA)*100, 1)}% of total cells, {np.round((nVOTselectiveTeA/nSoundResponsiveTeA)*100, 1)}% of sound responsive, {np.round((nVOTselectiveTeA/nSpeechResponsiveTeA)*100, 1)}% of speech responsive), n FT selective = {nFTselectiveTeA}, ({np.round((nFTselectiveTeA/nTotalTeA)*100, 1)}% of total cells, {np.round((nFTselectiveTeA/nSoundResponsiveTeA)*100, 1)}% of sound responsive, {np.round((nFTselectiveTeA/nSpeechResponsiveTeA)*100, 1)}% of speech responsive)')
    print('--Frac VOT selective of total cells --')
    print(f'AudP vs AudD p = {np.round(pvalFracVotSelective_AudPvsAudD,3)}')
    print(f'AudP vs AudV p = {np.round(pvalFracVotSelective_AudPvsAudV,3)}')
    print(f'AudP vs TeA p = {np.round(pvalFracVotSelective_AudPvsTeA,3)}')
    print(f'AudD vs AudV p = {np.round(pvalFracVotSelective_AudDvsAudV,3)}')
    print(f'AudD vs TeA p = {np.round(pvalFracVotSelective_AudDvsTeA,3)}')
    print(f'AudV vs TeA p = {np.round(pvalFracVotSelective_AudVvsTeA,3)}')
    print('--Frac VOT selective of sound responsive --')
    print(f'AudP vs AudD p = {np.round(pvalFracVotSelective_soundResponsiveAudPvsAudD,3)}')
    print(f'AudP vs AudV p = {np.round(pvalFracVotSelective_soundResponsiveAudPvsAudV,3)}')
    print(f'AudP vs TeA p = {np.round(pvalFracVotSelective_soundResponsiveAudPvsTeA,3)}')
    print(f'AudD vs AudV p = {np.round(pvalFracVotSelective_soundResponsiveAudDvsAudV,3)}')
    print(f'AudD vs TeA p = {np.round(pvalFracVotSelective_soundResponsiveAudDvsTeA,3)}')
    print(f'AudV vs TeA p = {np.round(pvalFracVotSelective_soundResponsiveAudVvsTeA,3)}')
    print('--Frac FT selective of total cells --')
    print(f'AudP vs AudD p = {np.round(pvalFracFtSelective_AudPvsAudD,3)}')
    print(f'AudP vs AudV p = {np.round(pvalFracFtSelective_AudPvsAudV,3)}')
    print(f'AudP vs TeA p = {np.round(pvalFracFtSelective_AudPvsTeA,3)}')
    print(f'AudD vs AudV p = {np.round(pvalFracFtSelective_AudDvsAudV,3)}')
    print(f'AudD vs TeA p = {np.round(pvalFracFtSelective_AudDvsTeA,3)}')
    print(f'AudV vs TeA p = {np.round(pvalFracFtSelective_AudVvsTeA,3)}')
    print('--Frac FT selective of sound responsive --')
    print(f'AudP vs AudD p = {np.round(pvalFracFtSelective_soundResponsiveAudPvsAudD,3)}')
    print(f'AudP vs AudV p = {np.round(pvalFracFtSelective_soundResponsiveAudPvsAudV,3)}')
    print(f'AudP vs TeA p = {np.round(pvalFracFtSelective_soundResponsiveAudPvsTeA,3)}')
    print(f'AudD vs AudV p = {np.round(pvalFracFtSelective_soundResponsiveAudDvsAudV,3)}')
    print(f'AudD vs TeA p = {np.round(pvalFracFtSelective_soundResponsiveAudDvsTeA,3)}')
    print(f'AudV vs TeA p = {np.round(pvalFracFtSelective_soundResponsiveAudVvsTeA,3)}')
    print('--quadrant-defined areas--')
    print(f'DP Quadrant: nTotal = {quadrantTotals[0]}, n Sound responsive = {np.sum(quadrantsSoundResponsive[0])}, n Speech responsive = {np.sum(quadrantsSpeechResponsive[0])}, n VOT selective = {np.sum(quadrantsVotSelective[0])} ({np.round((np.sum(quadrantsVotSelective[0])/quadrantTotals[0])*100,1)}% total, {np.round((np.sum(quadrantsVotSelective[0])/np.sum(quadrantsSoundResponsive[0]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsVotSelective[0])/np.sum(quadrantsSpeechResponsive[0]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[0])} ({np.round((np.sum(quadrantsFtSelective[0])/quadrantTotals[0])*100,1)}% total, {np.round((np.sum(quadrantsFtSelective[0])/np.sum(quadrantsSoundResponsive[0]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsFtSelective[0])/np.sum(quadrantsSpeechResponsive[0]))*100,1)}% of speech responsive)')
    print(f'DA Quadrant: nTotal = {quadrantTotals[1]}, n Sound responsive = {np.sum(quadrantsSoundResponsive[1])}, n Speech responsive = {np.sum(quadrantsSpeechResponsive[1])}, n VOT selective = {np.sum(quadrantsVotSelective[1])} ({np.round((np.sum(quadrantsVotSelective[1])/quadrantTotals[1])*100,1)}% total, {np.round((np.sum(quadrantsVotSelective[1])/np.sum(quadrantsSoundResponsive[1]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsVotSelective[1])/np.sum(quadrantsSpeechResponsive[1]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[1])} ({np.round((np.sum(quadrantsFtSelective[1])/quadrantTotals[1])*100,1)}% total, {np.round((np.sum(quadrantsFtSelective[1])/np.sum(quadrantsSoundResponsive[1]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsFtSelective[1])/np.sum(quadrantsSpeechResponsive[1]))*100,1)}% of speech responsive)')
    print(f'VP Quadrant: nTotal = {quadrantTotals[2]}, n Sound responsive = {np.sum(quadrantsSoundResponsive[2])}, n Speech responsive = {np.sum(quadrantsSpeechResponsive[2])}, n VOT selective = {np.sum(quadrantsVotSelective[2])} ({np.round((np.sum(quadrantsVotSelective[2])/quadrantTotals[2])*100,1)}% total, {np.round((np.sum(quadrantsVotSelective[2])/np.sum(quadrantsSoundResponsive[2]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsVotSelective[2])/np.sum(quadrantsSpeechResponsive[2]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[2])} ({np.round((np.sum(quadrantsFtSelective[2])/quadrantTotals[2])*100,1)}% total, {np.round((np.sum(quadrantsFtSelective[2])/np.sum(quadrantsSoundResponsive[2]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsFtSelective[2])/np.sum(quadrantsSpeechResponsive[2]))*100,1)}% of speech responsive)')
    print(f'VA Quadrant: nTotal = {quadrantTotals[3]}, n Sound responsive = {np.sum(quadrantsSoundResponsive[3])}, n Speech responsive = {np.sum(quadrantsSpeechResponsive[3])}, n VOT selective = {np.sum(quadrantsVotSelective[3])} ({np.round((np.sum(quadrantsVotSelective[3])/quadrantTotals[3])*100,1)}% total, {np.round((np.sum(quadrantsVotSelective[3])/np.sum(quadrantsSoundResponsive[3]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsVotSelective[3])/np.sum(quadrantsSpeechResponsive[3]))*100,1)}% of speech responsive), n FT selective = {np.sum(quadrantsFtSelective[3])} ({np.round((np.sum(quadrantsFtSelective[3])/quadrantTotals[3])*100,1)}% total, {np.round((np.sum(quadrantsFtSelective[3])/np.sum(quadrantsSoundResponsive[3]))*100,1)}% sound responsive, {np.round((np.sum(quadrantsFtSelective[3])/np.sum(quadrantsSpeechResponsive[3]))*100,1)}% of speech responsive)')
    print('--Frac VOT selective of total cells --')
    print(f'DP vs DA p = {np.round(quadrantComparFracVotSelective_allcells[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracVotSelective_allcells[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracVotSelective_allcells[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracVotSelective_allcells[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracVotSelective_allcells[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracVotSelective_allcells[5],3)}')
    print('--Frac VOT selective of sound responsive --')
    print(f'DP vs DA p = {np.round(quadrantComparFracVotSelective_soundResponsive[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracVotSelective_soundResponsive[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracVotSelective_soundResponsive[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracVotSelective_soundResponsive[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracVotSelective_soundResponsive[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracVotSelective_soundResponsive[5],3)}')
    print('--Frac FT selective of total cells --')
    print(f'DP vs DA p = {np.round(quadrantComparFracFtSelective_allcells[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracFtSelective_allcells[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracFtSelective_allcells[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracFtSelective_allcells[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracFtSelective_allcells[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracFtSelective_allcells[5],3)}')
    print('--Frac FT selective of sound responsive --')
    print(f'DP vs DA p = {np.round(quadrantComparFracFtSelective_soundResponsive[0],3)}')
    print(f'DP vs VP p = {np.round(quadrantComparFracFtSelective_soundResponsive[1],3)}')
    print(f'DP vs VA p = {np.round(quadrantComparFracFtSelective_soundResponsive[2],3)}')
    print(f'DA vs VP p = {np.round(quadrantComparFracFtSelective_soundResponsive[3],3)}')
    print(f'DA vs VA p = {np.round(quadrantComparFracFtSelective_soundResponsive[4],3)}')
    print(f'VP vs VA p = {np.round(quadrantComparFracFtSelective_soundResponsive[5],3)}')
    '''
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



if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
