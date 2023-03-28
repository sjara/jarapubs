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
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_cellLocations' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [6, 5] # In inches


fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
colorSpeechResp = cp.TangoPalette['ScarletRed2']
colorNotAud = cp.TangoPalette['Aluminium3']
colorSoundResp = cp.TangoPalette['Chocolate1']

# Load data
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath)

x_coords = figData['x_coord']
y_coords = figData['y_coord']
z_coords = figData['z_coord']
x_coords_jittered = figData['x_coords_jittered']
z_coords_jittered = figData['z_coords_jittered']
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
excludeCells = figData['excludeCells']

APtickLocs = np.array([156 ,176, 196, 216])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([200, 150, 100, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)
MLtickLocs = np.array([50, 150, 250, 350, 450])
MLtickLabels = np.round((250 - MLtickLocs)*0.025,1)

plt.figure()
gsMain = gridspec.GridSpec(1,1)
gsMain.update(left=0.1, right=0.8, top=0.8, bottom=0.1, wspace=0.4, hspace=0.3)

axDV_AP = plt.subplot(gsMain[0,0])
#axDV_ML = plt.subplot(gsMain[0,1])

plt.sca(axDV_AP)
nonResp = plt.scatter(z_coords_jittered[~speechResponsive | excludeCells], y_coords[~speechResponsive | excludeCells], c = colorNotAud, s=6)
soundResp = plt.scatter(z_coords_jittered[soundResponsive & ~speechResponsive], y_coords[soundResponsive & ~speechResponsive], c = colorSoundResp, s = 6)
speechResp = plt.scatter(z_coords_jittered[speechResponsive & ~excludeCells] , y_coords[speechResponsive & ~excludeCells], c = colorSpeechResp, s=6)
plt.xlim(155, 230)
plt.xticks(APtickLocs, APtickLabels)
plt.ylim(230,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Anterior - Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('Dorsal - Ventral (mm)', fontsize = fontSizeLabels)
plt.legend([nonResp, soundResp, speechResp], ['All cells', 'Sound responsive','Speech responsive'], loc = 'upper right', markerscale = 2 , bbox_to_anchor = (0.45,1.23))
axDV_AP.spines["right"].set_visible(False)
axDV_AP.spines["top"].set_visible(False)

'''
plt.sca(axDV_ML)
plt.scatter(x_coords_jittered[~speechResponsive | excludeCells], y_coords[~speechResponsive | excludeCells], c = colorNotAud, s=6)
plt.scatter(x_coords_jittered[speechResponsive & ~excludeCells], y_coords[speechResponsive & ~excludeCells], facecolor = colorSpeechResp, s=6)
plt.xlim(15,460)
plt.xticks(MLtickLocs, MLtickLabels)
plt.ylim(230,50)
plt.yticks(DVtickLocs, DVtickLabels)
plt.ylabel('Depth (mm)', fontsize = fontSizeLabels)
plt.xlabel('Lateral (mm)', fontsize = fontSizeLabels)
'''

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
