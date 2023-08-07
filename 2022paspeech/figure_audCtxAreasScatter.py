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
STATSUMMARY = 0
outputDir = settings.TEMP_OUTPUT_PATH
figFilename = 'audCtxAreas_scatter' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [3.5, 5] # In inches

figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

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
isCortical = figData['isCortical']
recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
APtickLocs = np.array([ 156 ,176, 196, 216, 236])
APtickLabels = np.round(-0.94 - (280-APtickLocs)*0.025,1)
DVtickLocs = np.array([210, 190, 170, 150, 130, 110, 90, 70, 50])
DVtickLabels = np.round((DVtickLocs-10)*0.025,1)


fig = plt.figure()
ax = fig.add_subplot()

plt.scatter(z_coords_jittered[recordingAreaName == audCtxAreas[0]], y_coords[recordingAreaName==audCtxAreas[0]], s=6)
plt.scatter(z_coords_jittered[recordingAreaName == audCtxAreas[1]], y_coords[recordingAreaName==audCtxAreas[1]], s=6)
plt.scatter(z_coords_jittered[recordingAreaName == audCtxAreas[2]], y_coords[recordingAreaName==audCtxAreas[2]], s=6)
plt.scatter(z_coords_jittered[recordingAreaName == audCtxAreas[3]], y_coords[recordingAreaName==audCtxAreas[3]], s=6)
plt.xlim(146, 246)
plt.xticks(APtickLocs, APtickLabels)
plt.ylim(220,40)
plt.yticks(DVtickLocs, DVtickLabels)
plt.xlabel('Posterior (mm)', fontsize = fontSizeLabels)
plt.ylabel('Ventral (mm)', fontsize = fontSizeLabels)
plt.legend(['AudP', 'AudD', 'AudV', 'TeA'])
ax.spines["right"].set_visible(False)
ax.spines["top"].set_visible(False)
ax.set_aspect('equal')


extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
