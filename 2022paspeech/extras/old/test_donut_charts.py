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
figFilename = 'figure_donutplots_responsive' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
figSize = [4, 4] # In inches

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizeTitles = figparams.fontSizeTitles
colorSpeechResp = cp.TangoPalette['ScarletRed2']
colorNotAud = cp.TangoPalette['Aluminium3']
colorSoundResp = cp.TangoPalette['Chocolate1']

# Load data
figDataFullPath = os.path.join(figDataDir,figDataFile)
figData = np.load(figDataFullPath, allow_pickle = True)

recordingAreaName = figData['recordingAreaName']
audCtxAreas = figData['audCtxAreas']
speechResponsive = figData['speechResponsive']
soundResponsive = figData['soundResponsive']
excludeCells = figData['excludeCells']

NspeechResponsiveByArea = np.empty(len(audCtxAreas))
NsoundResponsiveByArea = np.empty(len(audCtxAreas))
NtotalByArea = np.empty(len(audCtxAreas))

for indArea, thisArea in enumerate(audCtxAreas):
    NspeechResponsiveByArea[indArea] = np.sum(speechResponsive[(recordingAreaName == thisArea) & ~excludeCells])
    NsoundResponsiveByArea[indArea] = np.sum(soundResponsive[(recordingAreaName == thisArea) & ~excludeCells])
    NtotalByArea[indArea] = np.sum((recordingAreaName == thisArea) & ~excludeCells)

plt.figure()
gsMain = gridspec.GridSpec(2, 2)
axAudP = plt.subplot(gsMain[0,0])
axAudD = plt.subplot(gsMain[0,1])
axAudV = plt.subplot(gsMain[1,0])
axTeA = plt.subplot(gsMain[1,1])

circle1 = plt.Circle((0,0), 0.7, color = 'white')
circle2 = plt.Circle((0,0), 0.7, color = 'white')
circle3 = plt.Circle((0,0), 0.7, color = 'white')
circle4 = plt.Circle((0,0), 0.7, color = 'white')


plt.sca(axAudP)
plt.pie([NtotalByArea[0] - NsoundResponsiveByArea[0],  NsoundResponsiveByArea[0]-NspeechResponsiveByArea[0], NspeechResponsiveByArea[0]], colors = [colorNotAud, colorSoundResp, colorSpeechResp])
axAudP.add_artist(circle1)
plt.title(f'AudP, n = {int(NtotalByArea[0])}')

plt.sca(axAudD)
plt.pie([NtotalByArea[1] - NsoundResponsiveByArea[1],  NsoundResponsiveByArea[1]-NspeechResponsiveByArea[1], NspeechResponsiveByArea[1]], colors = [colorNotAud, colorSoundResp, colorSpeechResp])
axAudD.add_artist(circle2)
plt.title(f'AudD, n = {int(NtotalByArea[1])}')


plt.sca(axAudV)
plt.pie([NtotalByArea[2] - NsoundResponsiveByArea[2],  NsoundResponsiveByArea[2]-NspeechResponsiveByArea[2], NspeechResponsiveByArea[2]], colors = [colorNotAud, colorSoundResp, colorSpeechResp])
axAudV.add_artist(circle3)
plt.title(f'AudV, n = {int(NtotalByArea[2])}')


plt.sca(axTeA)
plt.pie([NtotalByArea[3] - NsoundResponsiveByArea[3],  NsoundResponsiveByArea[3]-NspeechResponsiveByArea[3], NspeechResponsiveByArea[3]], colors = [colorNotAud, colorSoundResp, colorSpeechResp])
axTeA.add_artist(circle4)
plt.title(f'TeA, n = {int(NtotalByArea[3])}')
plt.legend(['not sound responsive', 'sound responsive', 'speech responsive'], loc = 'lower right', bbox_to_anchor = (1.3, -0.35))

plt.show()


if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
