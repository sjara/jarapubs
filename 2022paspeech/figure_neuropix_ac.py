'''
This creates figure 2 for 2022paspeech:
 A. Cartoon of headfixed, awake mouse ephys
 B. Diagram of sound matrix
 C. Histology image of recording track
 D. Scatter plot of recording location of each cell. Add brain image in inkscape.
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
import figparams
import studyparams
from importlib import reload
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

reload(figparams)

FIGNAME = 'figure_neuropix_ac'
SAVE_FIGURE = 0
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_neuropix_ac' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [7, 4] # In inches

PANELS = [3,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
colorAudD = figparams.colors['audD']
colorAudP = figparams.colors['audP']
colorAudV = figparams.colors['audV']


labelPosX = [0.07, 0.3, 0.54] # Horiz position for panel labels
labelPosY = [0.93, 0.48]    # Vert position for panel labels

plt.figure()
gsMain = gridspec.GridSpec(1, 2, width_ratios=[0.3, 0.7])
gsMain.update(left=0.1, right=0.95, top=0.9, bottom=0.07, wspace=0.3, hspace=0.3)
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

y_coords = figData['y_coord'][~figData['excludeCells']]
z_coords = figData['z_coord'][~figData['excludeCells']]
recordingAreaName = figData['recordingAreaName'][~figData['excludeCells']]
isAudArea = figData['isAudArea'][~figData['excludeCells']]
vot_colors = figData['bestSelectivityIndexVot'][~figData['excludeCells']]
ft_colors = figData['bestSelectivityIndexFt'][~figData['excludeCells']]

## For now, plot on own figure

#plt.sca(axCellLocs)
'''
plt.scatter(z_coords[areaName == figData['audCtxAreas'][0]], y_coords[areaName == figData['audCtxAreas'][0]], c = colorAudP)
plt.scatter(z_coords[areaName == figData['audCtxAreas'][1]], y_coords[areaName == figData['audCtxAreas'][1]], c = colorAudD)
plt.scatter(z_coords[areaName == figData['audCtxAreas'][2]], y_coords[areaName == figData['audCtxAreas'][2]], c = colorAudV)
plt.legend(labels = ['AudP', 'AudD', 'AudV'], loc = 'lower left', bbox_to_anchor=(-0.3,0.8))
'''
#colors = plt.cm.viridis(vot_colors)
gsSelectivity = gridspec.GridSpec(1, 2)
axVot = plt.subplot(gsSelectivity[0,0])

colors = plt.cm.YlGn(vot_colors)
plt.scatter(z_coords, y_coords, c = colors)
plt.ylabel('D-V')
plt.xlabel('P-A')
axVOT.invert_yaxis()
plt.title('VOT selectivity by location')

axFt = plt.subplot(gsSelectivity[0,1])
colors = plt.cm.YlGn(ft_colors)
plt.scatter(z_coords, y_coords, c = colors)
plt.ylabel('D-V')
plt.xlabel('P-A')
axVOT.invert_yaxis()
plt.title('FT selectivity by location')
cbar = plt.colorbar(cmap = 'YlGn')
plt.show()

#axCellLocs.invert_yaxis()
#cbar = plt.colorbar(cmap = 'viridis')

#axFt = plt.subplot(gsSelectivity[0,1])



plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
