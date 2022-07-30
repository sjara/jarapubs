"""
Figure showing recording sites.
"""

import sys
sys.path.append('..')
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
outputDir = '/tmp/'
figFilename = 'fig2_locations' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
#figSize = [3.35, 2.4] # In inches
figSize = [7, 2.4] # In inches

'''
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel
'''

figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'sj_am_tuning_20220129.h5')
celldb = celldatabase.load_hdf(dbPath)

N_FREQ = 16 # HARDCODED
N_RATE = 11 # HARDCODED

# -- Find sound responsive cells --
toneResponsive = celldb['toneMinPval'] < (0.05/N_FREQ)
amOnsetResponsive = celldb['amMinPvalOnset'] < (0.05/N_RATE)
amSustainResponsive = celldb['amMinPvalSustain'] < (0.05/N_RATE)
soundResponsive = toneResponsive | amOnsetResponsive | amSustainResponsive

laserBase = celldb.laserBaseRate200
laserResp = celldb.laserRespRate50

# -- Exclude laser sessions with very low firing rate --
if 1:
    laserThFR = 1.0 #0.5 # Threshold in spk/s
    highFRlaserSession = (laserBase > laserThFR) | (laserResp > laserThFR)
    #highFRlaserSession = (laserBase > laserThFR)
    print(f'\nWARNING: Excluding cell with laser-session firing rate < {laserThFR} spk/s')
else:
    highFRlaserSession = True

cellsWithTone = ~np.isnan(celldb['toneMinPval'])
cellsWithAM = ~np.isnan(celldb['amMinPvalOnset'])
cellsWithSoundSessions = cellsWithTone | cellsWithAM
indD1 = ((celldb.laserPvalB200R50 < 0.01) & (laserResp>laserBase)) & highFRlaserSession
indND1 = ((celldb.laserPvalB200R50 > 0.01) | (laserResp<=laserBase)) & highFRlaserSession

# -- Select only high enough firing rate --
if 0:
    thFR = 1.0 # Threshold in spk/s
    highFR = ( (celldb.toneFiringRateBaseline > thFR) | (celldb.toneFiringRateBest > thFR) |
               (celldb.amFiringRateBaseline > thFR) | (celldb.amFiringRateBestOnset > thFR) |
               (celldb.amFiringRateBestSustain > thFR))
    indD1 = indD1 & highFR
    indND1 = indND1 & highFR
    print(f'\nWARNING: From the total, we analyze only cells with firing rate > {thFR} spk/s\n')

# -- Select sound-responsive --
respTypes = ['responsive','unresponsive']
indD1resp = indD1 & soundResponsive
indND1resp = indND1 & soundResponsive
indD1unresp = indD1 & ~soundResponsive
indND1unresp = indND1 & ~soundResponsive

# celldb.recordingSiteName.unique()
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

celldbCopy = celldb.copy()

#selectedSlices = np.array([244, 260, 270])
selectedSlices = np.array([266, 250])
zCoords = celldb.z_coord.to_numpy()
closestCoord = np.argmin(np.abs(zCoords[:,np.newaxis]-selectedSlices), axis=1)
newZcoord = selectedSlices[closestCoord].astype('float')
newZcoord[np.isnan(zCoords)] = np.nan
celldbCopy['z_coord'] = newZcoord

#celldbD1 = celldbCopy[indD1 & ~excludedByArea]
#celldbND1 = celldbCopy[indND1 & ~excludedByArea]
#celldbAll = celldbCopy[(indD1|indND1) & ~excludedByArea]

cellTypes = ['D1', 'ND1']
cellTypeLabel = ['D1', 'non-D1']
#indsEachType = [indD1, indND1]
indsEachType = [[indD1resp, indND1resp], [indD1unresp, indND1unresp]]

#colors = {'D1':cp.TangoPalette['SkyBlue1'],
#          'ND1':cp.TangoPalette['ScarletRed1']}
colors = {'D1':'#005fd8',
          'ND1':cp.TangoPalette['ScarletRed1']}

fontSizePanel = figparams.fontSizePanel
fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
panelLabels = [['A', 'B'],['C', 'D']]
labelPosX = [[0.01, 0.25], [0.54, 0.78]]   # Horiz position for panel labels
labelPosY = 0.93    # Vert position for panel labels
respLabelPosx = [[0.1125, 0.355], [0.645, 0.89]]

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

#gsMain = gridspec.GridSpec(2, 4)
gsMain = gridspec.GridSpec(1, 2)
gsMain.update(left=-0.01, right=1.01, top=1.0, bottom=0.05, wspace=0.1, hspace=0.0)

for indResp, responsiveness in enumerate(respTypes):
    gsResp = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=gsMain[indResp], wspace=0.0, hspace=0.0,)
    for indType, cellType in enumerate(cellTypes):
        siteColor = figparams.colors[cellType]
        #siteColor = colors[cellType]
        for inds, oneSlice in enumerate(selectedSlices):
            #thisAx = plt.subplot(gsMain[inds, indType + 2*indResp])
            thisAx = plt.subplot(gsResp[inds, indType])
            thisImgFile = os.path.join(figuresDataDir, f'atlas{oneSlice}_whiteBG.png')
            img = plt.imread(thisImgFile)
            plt.imshow(img, cmap='gray', vmax=255)
            theseCells = (celldbCopy.z_coord==oneSlice) & ~excludedByArea & indsEachType[indResp][indType]
            xCoords = celldb.x_coord[theseCells]
            yCoords = celldb.y_coord[theseCells]
            #plt.plot(xCoords, yCoords, '.', color=siteColor, ms=2)
            plt.plot(xCoords, yCoords, 'o', mec=siteColor, mfc='none', mew=0.5, ms=2)
            plt.axis(False)
        thisAx.annotate(panelLabels[indResp][indType], xy=(labelPosX[indResp][indType],labelPosY),
                        xycoords='figure fraction',
                        fontsize=fontSizePanel, fontweight='bold')
        respLabel = f'{cellTypeLabel[indType]} {respTypes[indResp]}'
        thisAx.annotate(respLabel, xy=(respLabelPosx[indResp][indType], 0.011),
                        xycoords='figure fraction', ha='center',
                        fontsize=fontSizeTicks, fontweight='normal')

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
