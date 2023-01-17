"""
Figure showing recording sites.
"""

import sys
import studyparams
import figparams
import os
import numpy as np
from jaratoolbox import settings
from jaratoolbox import celldatabase
from jaratoolbox import extraplots
from jaratoolbox import histologyanalysis as ha
from jaratoolbox import colorpalette as cp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from importlib import reload
reload(extraplots)


SAVE_FIGURE = 1
outputDir = 'C:\\Users\\jenny\\tmp'
figFilename = 'cell_locations' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
#figSize = [3.35, 2.4] # In inches
figSize = [7, 2.4] # In inches

'''
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
'''

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, figFilename)
databaseDir = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(databaseDir, 'fulldb_paspeech_speech_tuning_allcells.h5')
celldb = celldatabase.load_hdf(dbPath)

audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)

N_FREQ = 16 # HARDCODED
N_RATE = 11 # HARDCODED
N_SPEECH = 12

# -- Find sound responsive cells --
toneResponsive = celldb['toneMinPval'] < (0.05/N_FREQ)
amOnsetResponsive = celldb['amMinPvalOnset'] < (0.05/N_RATE)
amSustainResponsive = celldb['amMinPvalSustain'] < (0.05/N_RATE)
speechResponsive = celldb['speechMinPvalOnset'] < (0.05/N_SPEECH)
soundResponsive = toneResponsive | amOnsetResponsive | amSustainResponsive | speechResponsive


# -- Select only high enough firing rate --
if 0:
    exclusionCriterion = 4  # sp/s
    maxFiringRateSpeechEvoked = np.max([celldb.maxFiringRate_FT_VOTmin, celldb.maxFiringRate_FT_VOTmax, celldb.maxFiringRate_VOT_FTmin, celldb.maxFiringRate_VOT_FTmax], 0)
    exclude_evoked = maxFiringRateSpeechEvoked < exclusionCriterion
    exclude_baseline = celldb.speechFiringRateBaseline < exclusionCriterion
    excludeCells = exclude_evoked & exclude_baseline
    print(f'\nWARNING: From the total, we analyze only cells with firing rate > {exclusionCriterion} spk/s\n')

# -- Select sound-responsive --

respTypes = ['responsive','unresponsive']
'''
indD1resp = indD1 & soundResponsive
indND1resp = indND1 & soundResponsive
indD1unresp = indD1 & ~soundResponsive
indND1unresp = indND1 & ~soundResponsive
'''
indResp = soundResponsive #& ~excludeCells
indUnresp = ~soundResponsive #& excludeCells
# celldb.recordingSiteName.unique()
'''
areasToExclude = ['Supplemental somatosensory area, layer 6b',
                  'Supplemental somatosensory area, layer 6a',
                  'Visceral area, layer 6a',
                  'Primary somatosensory area, barrel field, layer 6b',
                  'Primary somatosensory area, barrel field, layer 6a',
                  'Primary somatosensory area, nose, layer 5']

excludedByArea = np.full(len(celldb), False, dtype='bool')
for inda, oneArea in enumerate(areasToExclude):
    celldb.recordingSiteName==oneArea
    excludedByArea[celldb.recordingSiteName==oneArea] = True
'''

areasToInclude = audCtxAreas
includeByArea = np.full(len(celldb), False, dtype = 'bool')
for indArea, thisArea in enumerate(audCtxAreas):
    celldb.recordingAreaName == thisArea
    includeByArea[celldb.recordingAreaName == thisArea] = True


celldbCopy = celldb.copy()

#selectedSlices = np.array([244, 260, 270])
selectedSlices = np.array([185, 200, 215])
zCoords = celldb.z_coord.to_numpy()
closestCoord = np.argmin(np.abs(zCoords[:,np.newaxis]-selectedSlices), axis=1)
newZcoord = selectedSlices[closestCoord].astype('float')
newZcoord[np.isnan(zCoords)] = np.nan
celldbCopy['z_coord'] = newZcoord

celldbAudP = celldbCopy[recordingAreaName == audCtxAreas[0]]
celldbAudD = celldbCopy[recordingAreaName == audCtxAreas[1]]
celldbAudV = celldbCopy[recordingAreaName == audCtxAreas[2]]
#celldbD1 = celldbCopy[indD1 & ~excludedByArea]
#celldbND1 = celldbCopy[indND1 & ~excludedByArea]
#celldbAll = celldbCopy[(indD1|indND1) & ~excludedByArea]

#cellTypes = ['D1', 'ND1']
#cellTypeLabel = ['D1', 'non-D1']
#indsEachType = [indD1, indND1]
#indsEachType = [[indD1resp, indND1resp], [indD1unresp, indND1unresp]]
indsEachType = [indResp, indUnresp]

#colors = {'D1':cp.TangoPalette['SkyBlue1'],
#          'ND1':cp.TangoPalette['ScarletRed1']}
#colors = {'D1':'#005fd8',
#          'ND1':cp.TangoPalette['ScarletRed1']}

fontSizePanel = figparams.fontSizePanel
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
panelLabels = [['A', 'B'],['C', 'D']]
labelPosX = [[0.01, 0.25], [0.54, 0.78]]   # Horiz position for panel labels
labelPosY = 0.93    # Vert position for panel labels
respLabelPosx = [[0.1125, 0.355], [0.645, 0.89]]

#aa = ha.AllenAverageCoronalAtlas()

# -- Plot results --
fig = plt.gcf()
fig.clf()
#fig.set_facecolor('w')

#gsMain = gridspec.GridSpec(2, 4)
gsMain = gridspec.GridSpec(1, 3)
gsMain.update(left=-0.01, right=1.01, top=1.0, bottom=0.05, wspace=0.1, hspace=0.0)

#for indResp, responsiveness in enumerate(respTypes):
#gsResp = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsMain, wspace=0.0, hspace=0.0,)
#for indType, cellType in enumerate(cellTypes):
#    siteColor = figparams.colors[cellType]
#    #siteColor = colors[cellType]
for inds, oneSlice in enumerate(selectedSlices):
    #thisAx = plt.subplot(gsMain[inds, indType + 2*indResp])
    thisAx = plt.subplot(gsMain[inds])
    thisImgFile = os.path.join(figuresDataDir, f'atlas{oneSlice}_blackBG.png')
    img = plt.imread(thisImgFile)
    #plt.imshow(img, cmap='gray', vmax=255)
    plt.imshow(img, cmap='gray') #, vmax=255)
    theseCells = (celldbCopy.z_coord==oneSlice) & includeByArea #& indsEachType[indResp]
    xCoords = celldb.x_coord[theseCells]
    yCoords = celldb.y_coord[theseCells]
    plt.plot(xCoords, yCoords, 'o', mfc='none', mew=0.5, ms=3)
    plt.axis(False)
    #respCells = (celldbCopy.z_coord==oneSlice) & includeByArea & indsEachType[indResp]
    #xCoords = celldb.x_coord[respCells]
    #yCoords = celldb.y_coord[respCells]
    #plt.plot(xCoords, yCoords, 'o', mec = cp.TangoPalette['ScarletRed1'], mfc='none', mew=0.5, ms=3)
    #plt.axis(False)


#thisAx.annotate(panelLabels[indResp], xy=(labelPosX[indResp],labelPosY),
#                xycoords='figure fraction',
#                fontsize=fontSizePanel, fontweight='bold')
#respLabel = f'{respTypes[indResp]}'
#thisAx.annotate(respLabel, xy=(respLabelPosx[indResp], 0.011),
#                xycoords='figure fraction', ha='center',
#                fontsize=fontSizeTicks, fontweight='normal')

plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir, dpi=300)


sys.exit()
# -- Use this to explore all sites --
if 0:
    aa = ha.AllenAverageCoronalAtlas()
    #aa.add_points_from_db(celldbAll)
    #aa.add_points_from_db(celldbAll)
    aa.add_points_from_db(celldbND1)
    aa.show_all_sites()


sys.exit()
# -- USE THIS CODE TO SAVE IMAGES --
inds = 266 #250
maxVal = 300  # Max pixel intensity (to scale to 8-bit)
plt.imsave(f'/tmp/atlas{inds}_blackBG.png', aa.atlas[:,:,inds], cmap='gray', vmax=maxVal)
'''
Then in GIMP:
FuzzySelect:  Radius=4, Threshold=10, FeatherEdges
(I created a white layer as background)
'''



#from PIL import Image
#imgArray = (255 * aa.atlas[:,:,inds]/maxVal).astype(int8)
#imgArray = aa.atlas[:,:,inds]
#img = Image.fromarray(imgArray)
#img.save(f'/tmp/atlas{inds}_blackgb.png')

''' '''
