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


SAVE_FIGURE = 0
AREA_OVERLAY = 0 #note, need to be in AllenSDK virtual env to include area overlay
#outputDir = 'C:\\Users\\jenny\\tmp'
outputDir = '/tmp/'
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

#selectedSlices = np.array([185, 200, 215])
selectedSlices = np.array([200])
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

fontSizePanel = figparams.fontSizePanel
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
panelLabels = [['A', 'B'],['C', 'D']]
labelPosX = [[0.01, 0.25], [0.54, 0.78]]   # Horiz position for panel labels
labelPosY = 0.93    # Vert position for panel labels
respLabelPosx = [[0.1125, 0.355], [0.645, 0.89]]


# -- Plot results --
if ~AREA_OVERLAY:
    fig = plt.gcf()
    fig.clf()
    #fig.set_facecolor('w')


    gsMain = gridspec.GridSpec(1, 3)
    gsMain.update(left=-0.01, right=1.01, top=1.0, bottom=0.05, wspace=0.1, hspace=0.0)

    for inds, oneSlice in enumerate(selectedSlices):
        thisAx = plt.subplot(gsMain[inds])
        thisImgFile = os.path.join(figuresDataDir, f'atlas{oneSlice}_blackBG.png')
        img = plt.imread(thisImgFile)
        #plt.imshow(img, cmap='gray', vmax=255)
        plt.imshow(img, cmap='gray')
        theseCells = (celldbCopy.z_coord==oneSlice) & includeByArea
        xCoords = celldb.x_coord[theseCells]
        yCoords = celldb.y_coord[theseCells]
        plt.plot(xCoords, yCoords, 'o', mfc='none', mew=0.5, ms=3)
        plt.axis(False)

    plt.show()

    if SAVE_FIGURE:
        extraplots.save_figure(figFilename, figFormat, figSize, outputDir, dpi=300)


#sys.exit()
# -- Use this to explore all sites --
if AREA_OVERLAY:
    aa = ha.AllenAverageCoronalAtlas()
    for inds, oneSlice in enumerate(selectedSlices):
        aa.get_slice(oneSlice)
        theseCells = (celldbCopy.z_coord==oneSlice) & includeByArea
        aa.add_points_from_db(celldbCopy[theseCells])
        aa.show_all_sites(areas=['AUDp', 'AUDv', 'AUDd', 'AUDpo'])

    if SAVE_FIGURE:
        extraplots.save_figure(figFilename, figFormat, figSize, outputDir, dpi=300)


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
