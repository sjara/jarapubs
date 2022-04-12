import os
import numpy as np
import nrrd
from matplotlib import pyplot as plt
from jaratoolbox import histologyanalysis as ha
from jaratoolbox import settings
reload(settings)
from jaratoolbox import extraplots
import pandas as pd
import re
import figparams

FIGNAME = 'figure_noise_laser'
figDir = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, FIGNAME)
# figDir = '/tmp'

dbPath = os.path.join(settings.FIGURES_DATA_PATH, figparams.STUDY_NAME, 'celldatabase_ALLCELLS.h5')
db = pd.read_hdf(dbPath, key='dataframe')

goodISI = db.query('isiViolations<0.02 or modifiedISI<0.02')
goodShape = goodISI.query('spikeShapeQuality > 2')
goodLaser = goodShape.query('autoTagged==1')

acCells = goodLaser.groupby('brainArea').get_group('rightAC')

atlasPath = os.path.join(settings.ALLEN_ATLAS_DIR, 'coronal_average_template_25.nrrd')
atlasData = nrrd.read(atlasPath)
atlas = atlasData[0]

#Hardcoded stuff
atlasWidth = 456
atlasHeight = 320

HISTOLOGY_PATH = settings.HISTOLOGY_PATH

plt.close('all')

#This splits the penetrations into 2 groups with similar atlas z values
zBound = 207

numFigs = 2
figs = []
axs = []
for indFig in range(numFigs):
    fig, ax = plt.subplots()
    ax.set_xlim([0, atlasWidth])
    ax.set_ylim([0, atlasHeight])
    ax.invert_yaxis()
    figs.append(fig)
    axs.append(ax)

lineColor = 'r'
lineWidth = 0.5
markerColor = 'r'
markerSize = 1

# fig1, ax1 = plt.subplots()

for indRow, dbRow in acCells.iterrows():

    # #For only 1 z section
    # tractAx = ax1
    # tractFig = fig1
    # sliceNum = 207

    if dbRow['cellZ']<zBound:
        tractAx = axs[0]
        tractFig = figs[0]
        sliceNum = 202
    elif (dbRow['cellZ']>zBound):
        tractAx = axs[1]
        tractFig = figs[1]
        sliceNum = 209

    # colors = {'pinp015':'r', 'pinp016':'c', 'pinp017':'g', 'pinp018':'b'}
    # markerColor = colors[dbRow['subject']]

    plt.sca(tractAx)
    plt.figure(tractFig.number)

    tractAx.imshow(np.rot90(atlas[:,:,sliceNum], -1), 'gray')
    plt.hold(1)
    # tractAx.plot([tract['tipCoords'][0], tract['brainSurfCoords'][0]], [tract['tipCoords'][1], tract['brainSurfCoords'][1]],
    #             color=lineColor, linewidth=lineWidth)
    # print tract['brainSurfCoords']
    # for coordPair in tract['siteCoords']:
    tractAx.plot(dbRow['cellX'], dbRow['cellY'], '.', mfc=markerColor, mec=markerColor, ms=markerSize)
    plt.hold(1)
    # plt.waitforbuttonpress()
    # print tract


# outputDir = '/mnt/jarahubdata/papers/2018thstr/figures/figure_noise_laser/'
outputDir = '/tmp'
for indFig, fig in enumerate(figs):
    plt.figure(fig.number)
    extraplots.save_figure('AC_tracts{}'.format(indFig+1), 'svg', (5, 3), outputDir=outputDir)
plt.show()

# tracts = [
# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1146., 'recordingDepths':[894, 1146],
#  'atlasZ':202, 'recordingTract':'medialDiI_shank1'},

# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1146., 'recordingDepths':[1146],
#  'atlasZ':205, 'recordingTract':'medialDiI_shank2'},

# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1146., 'recordingDepths':[957],
#  'atlasZ':209, 'recordingTract':'medialDiI_shank4'},

# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1503., 'recordingDepths':[957, 1175, 1275, 1378],
#  'atlasZ':201, 'recordingTract':'DiD_shank1'},

# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1503., 'recordingDepths':[975, 1087],
#  'atlasZ':207, 'recordingTract':'DiD_shank3'},

# {'subject':'pinp015', 'brainArea':'AC',
#  'maxDepth':1503., 'recordingDepths':[975, 1087, 1175],
#  'atlasZ':210, 'recordingTract':'DiD_shank4'},

# {'subject':'pinp016', 'brainArea':'AC',
#  'maxDepth':2051., 'recordingDepths':[2051],
#  'atlasZ':214, 'recordingTract':'lateralDiI_shank1'},

# {'subject':'pinp016', 'brainArea':'AC',
#  'maxDepth':2051., 'recordingDepths':[1904],
#  'atlasZ':202, 'recordingTract':'lateralDiI_shank3'},

# {'subject':'pinp016', 'brainArea':'AC',
#  'maxDepth':2051., 'recordingDepths':[1153, 1904, 2051],
#  'atlasZ':197, 'recordingTract':'lateralDiI_shank4'},

# {'subject':'pinp016', 'brainArea':'AC',
#  'maxDepth':2091., 'recordingDepths':[1143, 1338],
#  'atlasZ':193, 'recordingTract':'extraLateralDiD_shank3'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1338., 'recordingDepths':[1143],
#  'atlasZ':200, 'recordingTract':'medialDiI_shank1'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1338., 'recordingDepths':[1143, 1338],
#  'atlasZ':210, 'recordingTract':'medialDiI_shank2'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1338., 'recordingDepths':[1338],
#  'atlasZ':217, 'recordingTract':'medialDiI_shank3'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1338., 'recordingDepths':[1143, 1247, 1338],
#  'atlasZ':228, 'recordingTract':'medialDiI_shank4'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1604., 'recordingDepths':[1281, 1604],
#  'atlasZ':198, 'recordingTract':'medialDiD_shank1'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1604., 'recordingDepths':[1281, 1518],
#  'atlasZ':205, 'recordingTract':'medialDiD_shank2'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1604., 'recordingDepths':[1414, 1518],
#  'atlasZ':211,'recordingTract':'medialDiD_shank3'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1604., 'recordingDepths':[1281, 1518, 1604],
#  'atlasZ':218,'recordingTract':'medialDiD_shank4'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1525., 'recordingDepths':[1401],
#  'atlasZ':213,'recordingTract':'lateralDiI_shank1'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1525., 'recordingDepths':[1401, 1525],
#  'atlasZ':215,'recordingTract':'lateralDiI_shank3'},

# {'subject':'pinp017', 'brainArea':'AC',
#  'maxDepth':1525., 'recordingDepths':[1401],
#  'atlasZ':220,'recordingTract':'lateralDiI_shank4'},

# {'subject':'pinp018', 'brainArea':'AC',
#  'maxDepth':1136., 'recordingDepths':[1016],
#  'atlasZ':200,'recordingTract':'DiD_shank2'},

# {'subject':'pinp018', 'brainArea':'AC',
#  'maxDepth':1136., 'recordingDepths':[966],
#  'atlasZ':201,'recordingTract':'DiD_shank4'}
#     ]

# for tract in tracts:
#     registrationFolder = 'registration{}'.format(tract['brainArea'])
#     filenameSVG = os.path.join(HISTOLOGY_PATH, tract['subject'], registrationFolder, '{}.svg'.format(tract['recordingTract']))
#     tree = ETree.parse(filenameSVG)
#     root=tree.getroot()
#     paths = root.findall('{http://www.w3.org/2000/svg}path')
#     if len(paths)!=1:
#         raise ValueError('The SVG file must contain exactly 1 path')
#     pathCoords = paths[0].attrib['d']

#     # FIXME:This is the correct way to do it, but not working for me and I'm fed up
#     # reString = r'M (\d+\.*\d*),(\d+\.*\d*) (\d+\.*\d*),(\d+\.*\d*)'
#     # coordStrings = re.findall(reString, pathCoords, flags=re.IGNORECASE)

#     stackedCoordStrings = [coords.split(',') for coords in pathCoords.split(' ')[1:]]
#     coordStrings = []
#     for coordPair in stackedCoordStrings:
#         for coord in coordPair:
#             coordStrings.append(coord)

#     # if len(coordStrings)==0:
#     if len(coordStrings)!=4:
#         raise ValueError('The path does not have the correct format. You probably did not double click for this tract')
#     # tractCoords = coordStrings[0]
#     # tractCoords = map(float, tractCoords)
#     tractCoords = map(float, coordStrings)
#     # line = Line(tractCoords)

#     tipCoords = [tractCoords[0], tractCoords[1]]
#     brainSurfCoords = [tractCoords[2], tractCoords[3]]

#     #Was doing this because sometimes Inkscape will use relative path coords and it will break things. 
#     import ipdb
#     if any(np.array(tipCoords)<1):
#         ipdb.set_trace()
#     if any(np.array(brainSurfCoords)<1):
#         ipdb.set_trace()

#     if tipCoords[1] < brainSurfCoords[1]:
#         # raise ValueError('The brain surface is deeper than the tip!')
#         print "Bad coords, skipping {}".format(tract)
#         continue

#     siteFracFromSurface = np.array(tract['recordingDepths'])/float(tract['maxDepth'])
#     siteCoords = [tract_fraction(tipCoords, brainSurfCoords, fracFromSurface) for fracFromSurface in siteFracFromSurface]
#     tract.update({'siteCoords':siteCoords, 'tipCoords':tipCoords, 'brainSurfCoords':brainSurfCoords})

# for tract in tracts:

#     # #For only 1 z section
#     # tractAx = ax1
#     # tractFig = fig1
#     # sliceNum = 207

#     if tract['atlasZ']<zBound:
#         tractAx = axs[0]
#         tractFig = figs[0]
#         sliceNum = 202
#     elif (tract['atlasZ']>zBound):
#         tractAx = axs[1]
#         tractFig = figs[1]
#         sliceNum = 209

#     plt.sca(tractAx)
#     plt.figure(tractFig.number)

#     tractAx.imshow(np.rot90(atlas[:,:,sliceNum], -1), 'gray')
#     plt.hold(1)
#     # tractAx.plot([tract['tipCoords'][0], tract['brainSurfCoords'][0]], [tract['tipCoords'][1], tract['brainSurfCoords'][1]],
#     #             color=lineColor, linewidth=lineWidth)
#     print tract['brainSurfCoords']
#     for coordPair in tract['siteCoords']:
#         tractAx.plot(coordPair[0], coordPair[1], '.', mfc=markerColor, mec=markerColor, ms=1)
#     plt.hold(1)
#     # plt.waitforbuttonpress()
#     # print tract
