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
from importlib import reload
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.font_manager import findfont, FontProperties
font = findfont(FontProperties(family = ['Helvetica']))
reload(figparams)

FIGNAME = 'selectivityIndices'
SAVE_FIGURE = 0
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_neuropix_methods' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 5] # In inches

PANELS = [3,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
fontSizeTitles = figparams.fontSizeTitles
colorSpeechResp = cp.TangoPalette['ScarletRed2']
colorNotAud = cp.TangoPalette['Aluminium3']
colorSoundResp = cp.TangoPalette['SkyBlue2']


labelPosX = [0.07, 0.3, 0.53] # Horiz position for panel labels
labelPosY = [0.93, 0.48]    # Vert position for panel labels

plt.figure()
gsMain = gridspec.GridSpec(1, 2, width_ratios=[0.3, 0.7])
gsMain.update(left=0.1, right=0.93, top=0.9, bottom=0.1, wspace=0.3, hspace=0.3)
gsLeft = gsMain[0].subgridspec(2,1)

axCartoon = plt.subplot(gsLeft[0,0])
axCartoon.set_axis_off()

axSoundMatrix = plt.subplot(gsLeft[1,0])
axSoundMatrix.set_axis_off()

gsRight = gsMain[1].subgridspec(1,2, width_ratios = [0.3, 0.7])
axHist = plt.subplot(gsRight[0,0])
axHist.set_axis_off()

axCellLocs = plt.subplot(gsRight[0,1])


axCartoon.annotate('A', xy=(labelPosX[0],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axSoundMatrix.annotate('B', xy=(labelPosX[0],labelPosY[1]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axHist.annotate('C', xy=(labelPosX[1],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
axCellLocs.annotate('D', xy=(labelPosX[2],labelPosY[0]), xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')


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
excludeCells = figData['excludeCells']
recordingAreaName = figData['recordingAreaName']

APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)


plt.sca(axCellLocs)
nonResp = plt.scatter(z_coords_jittered[~(speechResponsive | excludeCells)], y_coords[~(speechResponsive | excludeCells)], c = colorNotAud, s=6)
soundResp = plt.scatter(z_coords_jittered[soundResponsive & ~speechResponsive], y_coords[soundResponsive & ~speechResponsive], c = colorSoundResp, s = 6)
speechResp = plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells] , y_coords[speechResponsive & ~excludeCells], c = colorSpeechResp, s=6)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Posterior from Bregma (mm)', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.legend([nonResp, soundResp, speechResp], ['All cells', 'Sound responsive','Speech responsive'], loc = 'upper center', markerscale = 2 , bbox_to_anchor = (0.5,1.1))
axCellLocs.spines["right"].set_visible(False)
axCellLocs.spines["top"].set_visible(False)


#plt.sca(axCellLocs)
'''
plt.scatter(z_coords[areaName == figData['audCtxAreas'][0]], y_coords[areaName == figData['audCtxAreas'][0]], c = colorAudP)
plt.scatter(z_coords[areaName == figData['audCtxAreas'][1]], y_coords[areaName == figData['audCtxAreas'][1]], c = colorAudD)
plt.scatter(z_coords[areaName == figData['audCtxAreas'][2]], y_coords[areaName == figData['audCtxAreas'][2]], c = colorAudV)
plt.legend(labels = ['AudP', 'AudD', 'AudV'], loc = 'lower left', bbox_to_anchor=(-0.3,0.8))
'''


#axCellLocs.invert_yaxis()
#cbar = plt.colorbar(cmap = 'viridis')

#axFt = plt.subplot(gsSelectivity[0,1])



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
