"""
Creates firing rate based selectivity indices for VOT and FT (maximum modulation by VOT for each FT) and compares these indices between cortical areas. Also checks following controls (comparing between areas): average of VOT SI between FTmin/FTmax, checking SI-VOT FTmin and SI-VOT FTmax separately, checking max(SI-VOT FTmin, SI-VOT FTmax) including all cells (not just cells that respond to speech sounds)
"""

import os
import sys
import numpy as np
from jaratoolbox import settings
from jaratoolbox import celldatabase
from scipy import stats
import studyparams
import studyutils

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
STATSUMMARY = 1
shuffledDataFile = 'data_shuffledSIs.npz'
boundariesDataFile = 'brain_areas_boundaries.npz'


if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)

figDataFullPath = os.path.join(figDataDir,figDataFile)

shuffledDataFullPath = os.path.join(figDataDir,shuffledDataFile)
shuffledData = np.load(shuffledDataFullPath, allow_pickle=True)
boundariesDataFullPath = os.path.join(figDataDir, boundariesDataFile)
boundariesData = np.load(boundariesDataFullPath)
#scriptFullPath = os.path.realpath(__file__)



databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']

audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area', 'Temporal association areas']
recordingAreaName = celldb.recordingAreaName.copy()
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)
isAudArea = celldb.recordingAreaName.str.contains('auditory area')
ctxAreas = ['Anterolateral visual area', 'Dorsal auditory area', 'Ectorhinal area', 'Endopiriform nucleus', 'Lateral visual area', 'Laterointermediate area', 'Perirhinal area', 'Posterior auditory area', 'Primary auditory area', 'Supplemental somatosensory area', 'Temporal association areas', 'Ventral auditory area', 'auditory radiation', 'optic radiation', 'corpus callosum', 'external capsule']

isCortical = np.zeros(len(isAudArea), dtype = bool)
for indArea, thisArea in enumerate(ctxAreas):
    isCortical[recordingAreaName == thisArea] = True


layersDeep = celldb.recordingSiteName.str.contains('layer 5|layer 6') & isCortical
layer4 =  celldb.recordingSiteName.str.contains('layer 4') & isCortical
layersSuperficial =  celldb.recordingSiteName.str.contains('layer 1|layer 2| layer 3') & isCortical
layersDeep[isCortical & ~(layersDeep|layer4|layersSuperficial)] = True



selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])



selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])


## -- exclude cells with low firing rates
exclusionCriterion = 5  # sp/s
maxFiringRateSpeechEvoked = np.max([celldb.maxFiringRate_FT_VOTmin, celldb.maxFiringRate_FT_VOTmax, celldb.maxFiringRate_VOT_FTmin, celldb.maxFiringRate_VOT_FTmax], 0)
excludeSpeech = (maxFiringRateSpeechEvoked < exclusionCriterion) & (celldb.speechFiringRateBaseline < exclusionCriterion)
excludeAM = (celldb.amFiringRateMaxOnset < exclusionCriterion) & (celldb.amFiringRateBaseline < exclusionCriterion)
excludeTone = (celldb.toneFiringRateBest < exclusionCriterion) & (celldb.amFiringRateBaseline < exclusionCriterion)

manualExcludedCells = [1534, 1532, 1147] #manually exclude cells that are accidental subsets of other cells due to Phy issues.
excludeSpeech[manualExcludedCells] = True
excludeAM[manualExcludedCells] = True
excludeTone[manualExcludedCells] = True


speechAlpha = 0.05/12
speechResponsive = np.array((celldb.speechMinPvalOnset < speechAlpha) & ~excludeSpeech)
amAlpha = 0.05/11 #11 AM frequencies presented
toneAlpha = 0.05/16 #16 pure tones presented
toneResponsive = np.array((celldb.toneMinPval < toneAlpha) & ~excludeTone)
amResponsive = np.array(((celldb.amMinPvalOnset < amAlpha) |
                    (celldb.amMinPvalSustain < amAlpha)) & ~excludeAM)
soundResponsive = np.array(toneResponsive|amResponsive|speechResponsive)


## -- select best index b/w conditions
ftCombos = np.array([selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax])
votCombos = np.array([selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax])


amOnsetSelective = np.array(celldb.amSelectivityPvalOnset < 0.05)
amSustainSelective = np.array(celldb.amSelectivityPvalSustain < 0.05)
amSelective = amOnsetSelective| amSustainSelective
toneSelective = np.array(celldb.toneSelectivityPval < 0.05)

amSelective[excludeAM] = np.nan
toneSelective[excludeTone] = np.nan

# -- Calculate pvalues for significant selectivity (permutation test)
#shuffledSI = shuffledData['shuffledSI']
#shuffledSiEachCell = np.concatenate((shuffledSI[0], shuffledSI[1], shuffledSI[2], shuffledSI[3], shuffledSI[4], shuffledSI[5], shuffledSI[6]))
#avgShuffledSiEachCell = np.mean(shuffledSiEachCell, 1)
#avgShuffledSI = np.nanmean(avgShuffledSiEachCell)

nCells = len(shuffledData['votSI_FTmax'])
votSI_FTmax = shuffledData['votSI_FTmax']
votSI_FTmin = shuffledData['votSI_FTmin']
shuffledSIVot_FTmax = shuffledData['shuffledSIVot_FTmax']
shuffledSIVot_FTmin = shuffledData['shuffledSIVot_FTmin']
ftSI_VOTmax = shuffledData['ftSI_VOTmax']
ftSI_VOTmin = shuffledData['ftSI_VOTmin']
shuffledSIFt_VOTmax = shuffledData['shuffledSIFt_VOTmax']
shuffledSIFt_VOTmin = shuffledData['shuffledSIFt_VOTmin']
maxRateVotFTmax = shuffledData['maxRateVotFTmax']
maxRateVotFTmin = shuffledData['maxRateVotFTmin']
maxRateFtVOTmax = shuffledData['maxRateFtVOTmax']
maxRateFtVOTmin = shuffledData['maxRateFtVOTmin']
minRateVotFTmax = shuffledData['minRateVotFTmax']
minRateVotFTmin = shuffledData['minRateVotFTmin']
minRateFtVOTmax = shuffledData['minRateFtVOTmax']
minRateFtVOTmin = shuffledData['minRateFtVOTmin']


pvalPermutationtestFt = np.ones(nCells)
pvalPermutationtestVot = np.ones(nCells)
bestSelectivityIndexFt = np.empty(nCells)
bestSelectivityIndexVot = np.empty(nCells)
shuffledVotBest = np.empty(np.shape(shuffledSIVot_FTmax))
shuffledFtBest = np.empty(np.shape(shuffledSIVot_FTmax))
whichFT = np.full(nCells, np.nan, dtype = int)
whichVOT = np.full(nCells, np.nan, dtype = int)

whichFT_labels = {0:'FT_VOTmin',1:'FT_VOTmax','FT_VOTmin':0,'FT_VOTmax':1}
whichVOT_labels = {0:'VOT_FTmin',1:'VOT_FTmax','VOT_FTmin':0,'VOT_FTmax':1}
frMinBestSI = 5

rateDiffBaselineVot_FTmax = np.abs(celldb.speechFiringRateBaseline - maxRateVotFTmax)
rateDiffBaselineVot_FTmin = np.abs(celldb.speechFiringRateBaseline - maxRateVotFTmin)
rateDiffBaselineFt_VOTmax = np.abs(celldb.speechFiringRateBaseline - maxRateFtVOTmax)
rateDiffBaselineFt_VOTmin = np.abs(celldb.speechFiringRateBaseline - maxRateFtVOTmin)


for indCell, thisCell in enumerate(votSI_FTmax):
    if (rateDiffBaselineVot_FTmax[indCell] > rateDiffBaselineVot_FTmin[indCell]) & (maxRateVotFTmax[indCell] >= frMinBestSI):
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmax[indCell,:] >=votSI_FTmax[indCell])
        bestSelectivityIndexVot[indCell] = votSI_FTmax[indCell]
        shuffledVotBest[indCell] = shuffledSIVot_FTmax[indCell]
        whichVOT[indCell] = 1
    elif (rateDiffBaselineVot_FTmin[indCell] > rateDiffBaselineVot_FTmax[indCell]) & (maxRateVotFTmin[indCell] >= frMinBestSI):
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmin[indCell,:] >=votSI_FTmin[indCell])
        bestSelectivityIndexVot[indCell] = votSI_FTmin[indCell]
        shuffledVotBest[indCell] = shuffledSIVot_FTmin[indCell]
        whichVOT[indCell] = 0
    elif (rateDiffBaselineVot_FTmax[indCell] > rateDiffBaselineVot_FTmin[indCell]) & (maxRateVotFTmax[indCell] < frMinBestSI) & (maxRateVotFTmin[indCell] >= frMinBestSI):
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmin[indCell,:] >=votSI_FTmin[indCell])
        bestSelectivityIndexVot[indCell] = votSI_FTmin[indCell]
        shuffledVotBest[indCell] = shuffledSIVot_FTmin[indCell]
        whichVOT[indCell] = 0
    elif (rateDiffBaselineVot_FTmin[indCell] > rateDiffBaselineVot_FTmax[indCell]) & (maxRateVotFTmin[indCell] < frMinBestSI) & (maxRateVotFTmax[indCell] >= frMinBestSI):
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmax[indCell,:] >=votSI_FTmax[indCell])
        bestSelectivityIndexVot[indCell] = votSI_FTmax[indCell]
        shuffledVotBest[indCell] = shuffledSIVot_FTmax[indCell]
        whichVOT[indCell] = 1
    else:
        bestSelectivityIndexVot[indCell] = votSI_FTmax[indCell]
        shuffledVotBest[indCell] = shuffledSIVot_FTmax[indCell]
        whichVOT[indCell] = 1

    if (rateDiffBaselineFt_VOTmax[indCell] > rateDiffBaselineFt_VOTmin[indCell]) & (maxRateFtVOTmax[indCell] >= frMinBestSI):
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmax[indCell,:] >=ftSI_VOTmax[indCell])
        bestSelectivityIndexFt[indCell] = ftSI_VOTmax[indCell]
        shuffledFtBest[indCell] = shuffledSIFt_VOTmax[indCell]
        whichFT[indCell] = 1
    elif (rateDiffBaselineFt_VOTmin[indCell] > rateDiffBaselineFt_VOTmax[indCell]) & (maxRateFtVOTmin[indCell] >= frMinBestSI):
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmin[indCell,:] >= ftSI_VOTmin[indCell])
        bestSelectivityIndexFt[indCell] = ftSI_VOTmin[indCell]
        shuffledFtBest[indCell] = shuffledSIFt_VOTmin[indCell]
        whichFT[indCell] = 0
    elif (rateDiffBaselineFt_VOTmax[indCell] > rateDiffBaselineFt_VOTmin[indCell]) & (maxRateFtVOTmax[indCell] < frMinBestSI) & (maxRateFtVOTmin[indCell] >= frMinBestSI):
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmin[indCell,:] >= ftSI_VOTmin[indCell])
        bestSelectivityIndexFt[indCell] = ftSI_VOTmin[indCell]
        shuffledFtBest[indCell] = shuffledSIFt_VOTmin[indCell]
        whichFT[indCell] = 0
    elif (rateDiffBaselineFt_VOTmin[indCell] > rateDiffBaselineFt_VOTmax[indCell]) & (maxRateFtVOTmin[indCell] < frMinBestSI) & (maxRateFtVOTmax[indCell] >= frMinBestSI):
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmax[indCell,:] >=ftSI_VOTmax[indCell])
        bestSelectivityIndexFt[indCell] = ftSI_VOTmax[indCell]
        shuffledFtBest[indCell] = shuffledSIFt_VOTmax[indCell]
        whichFT[indCell] = 1
    else:
        bestSelectivityIndexFt[indCell] = ftSI_VOTmax[indCell]
        shuffledFtBest[indCell] = shuffledSIFt_VOTmax[indCell]
        whichFT[indCell] = 1



## -- get number of selective cells for each feature
selectivityCriterion = 0.05
VOTselective = pvalPermutationtestVot < selectivityCriterion
FTselective = pvalPermutationtestFt < selectivityCriterion
mixedSelective = VOTselective & FTselective
singleSelective = np.logical_xor(VOTselective, FTselective)
notSelective = ~VOTselective & ~FTselective

y_coords = celldb.y_coord.copy()
z_coords = celldb.z_coord.copy()
x_coords = celldb.x_coord.copy()

z_coords_jittered = z_coords + np.random.randn(len(z_coords))
x_coords_jittered = x_coords + np.random.randn(len(x_coords))

bestFtSIbyArea = []
bestVotSIbyArea = []
soundResponsiveByArea = []
speechResponsiveByArea = []
VOTselectivebyArea = []
FTselectivebyArea = []
amSelectivebyArea = []
toneSelectivebyArea = []
excludeSpeechbyArea = []
mixedSelectivebyArea = []
singleSelectivebyArea = []
notSelectivebyArea = []
y_coordsbyArea = []
z_coordsbyArea = []

for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea.append(bestSelectivityIndexFt[recordingAreaName == thisArea])
    bestVotSIbyArea.append(bestSelectivityIndexVot[recordingAreaName == thisArea])
    soundResponsiveByArea.append(soundResponsive[recordingAreaName == thisArea])
    speechResponsiveByArea.append(speechResponsive[recordingAreaName == thisArea])
    amSelectivebyArea.append(amSelective[recordingAreaName == thisArea])
    toneSelectivebyArea.append(toneSelective[recordingAreaName == thisArea])
    excludeSpeechbyArea.append(excludeSpeech[recordingAreaName == thisArea])
    VOTselectivebyArea.append(VOTselective[recordingAreaName == thisArea])
    FTselectivebyArea.append(FTselective[recordingAreaName == thisArea])
    mixedSelectivebyArea.append(mixedSelective[recordingAreaName == thisArea])
    singleSelectivebyArea.append(singleSelective[recordingAreaName == thisArea])
    notSelectivebyArea.append(notSelective[recordingAreaName == thisArea])
    y_coordsbyArea.append(y_coords[recordingAreaName == thisArea])
    z_coordsbyArea.append(z_coords[recordingAreaName == thisArea])


## -- exclude low spike count cells
for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea[indArea] = bestFtSIbyArea[indArea][~excludeSpeechbyArea[indArea]]
    bestVotSIbyArea[indArea] = bestVotSIbyArea[indArea][~excludeSpeechbyArea[indArea]]
    soundResponsiveByArea[indArea] = soundResponsiveByArea[indArea][~excludeSpeechbyArea[indArea]]
    speechResponsiveByArea[indArea] = speechResponsiveByArea[indArea][~excludeSpeechbyArea[indArea]]
    amSelectivebyArea[indArea] = amSelectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    toneSelectivebyArea[indArea] = toneSelectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    VOTselectivebyArea[indArea] = VOTselectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    FTselectivebyArea[indArea] = FTselectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    mixedSelectivebyArea[indArea] = mixedSelectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    singleSelectivebyArea[indArea] = singleSelectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    notSelectivebyArea[indArea] = notSelectivebyArea[indArea][~excludeSpeechbyArea[indArea]]
    y_coordsbyArea[indArea] = y_coordsbyArea[indArea][~excludeSpeechbyArea[indArea]]
    z_coordsbyArea[indArea] = z_coordsbyArea[indArea][~excludeSpeechbyArea[indArea]]


# -- QUADRANTS
'''
audCtxAPbounds = np.array([-3.95,-1.9])
spreadAP = audCtxAPbounds[1]-audCtxAPbounds[0]
quadrantAPThreshold = audCtxAPbounds[0] + (spreadAP/2)
quadrantBoundsAP = np.array([audCtxAPbounds[0], quadrantAPThreshold, audCtxAPbounds[1]]) #in mm
quadrantBoundsAP_AApix = 280 + (quadrantBoundsAP + 0.94)/0.025 #in pixels
'''
extentAP = boundariesData['extentAP']
spreadAP = extentAP[1] - extentAP[0]
quadrantAPThreshold =  extentAP[0] + (spreadAP/2)
quadrantBoundsAP_AApix = np.array([extentAP[0], quadrantAPThreshold, extentAP[1]])
quadrantBoundsAP = studyutils.pix2mmAP(quadrantBoundsAP_AApix)

extentDV = boundariesData['extentDV']
spreadDV = extentDV[1] - extentDV[0]
quadrantDVThreshold =  extentDV[0] + (spreadDV/2)
quadrantBoundsDV_AApix = np.array([extentDV[0], quadrantDVThreshold, extentDV[1]])
quadrantBoundsDV = studyutils.pix2mmDV(quadrantBoundsDV_AApix)



quantilesDV = np.zeros([2, len(soundResponsive)], dtype = bool)
quantilesAP = np.zeros([2, len(soundResponsive)], dtype = bool)
for indBin, thisBin in enumerate(quantilesDV):
    thisQuantileDV = (y_coords >= quadrantBoundsDV_AApix[indBin]) & (y_coords < quadrantBoundsDV_AApix[indBin+1]) & isCortical & ~excludeSpeech
    thisQuantileAP = (z_coords >= quadrantBoundsAP_AApix[indBin]) & (z_coords < quadrantBoundsAP_AApix[indBin+1]) & isCortical & ~excludeSpeech
    quantilesAP[indBin] = thisQuantileAP
    quantilesDV[indBin] = thisQuantileDV

quadrantLabelsAP = ['posterior', 'anterior']
quadrantLabelsDV = ['dorsal', 'ventral']
quadrantLabels = ['dorsal posterior', 'dorsal anterior', 'ventral posterior', 'ventral anterior']

quadrantTotals = np.array(np.zeros(len(quadrantLabels), dtype=int))
quadrantsVOT = []
quadrantsFT = []
quadrantsVotSelective = []
quadrantsFtSelective = []
quadrantsSingleSelective = []
quadrantsMixedSelective = []
quadrantsSpeechResponsive = []
quadrantsSoundResponsive = []

a = -1
for indBinDV, thisQuantileDV in enumerate(quantilesDV):
    for indBinAP, thisQuantileAP in enumerate(quantilesAP):
        a = a+1
        quadrantTotals[a] = np.sum(thisQuantileDV & thisQuantileAP)
        quadrantsVOT.append(bestSelectivityIndexVot[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsFT.append(bestSelectivityIndexFt[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsVotSelective.append(VOTselective[thisQuantileAP & thisQuantileDV & speechResponsive])
        quadrantsFtSelective.append(FTselective[thisQuantileAP & thisQuantileDV & speechResponsive])
        quadrantsSingleSelective.append(singleSelective[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsMixedSelective.append(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive])
        quadrantsSpeechResponsive.append(speechResponsive[thisQuantileDV & thisQuantileAP])
        quadrantsSoundResponsive.append(soundResponsive[thisQuantileDV & thisQuantileAP])



# -- TEST EACH ANIMAL
CHECKBYANIMAL = 1
nAnimals = len(np.unique(celldb.subject))
quadrantTotalsByAnimal = np.empty([nAnimals, 4], dtype=int)
quadrantsSoundResponsiveByAnimal = np.empty([nAnimals, 4], dtype = int)
quadrantsSpeechResponsiveByAnimal = np.empty([nAnimals, 4], dtype = int)
quadrantsVotSelectiveByAnimal = np.empty([nAnimals, 4], dtype = int)
quadrantsFtSelectiveByAnimal = np.empty([nAnimals, 4], dtype = int)
quadrantsSingleSelectiveByAnimal = np.empty([nAnimals, 4], dtype = int)
quadrantsMixedSelectiveByAnimal = np.empty([nAnimals, 4], dtype = int)



if CHECKBYANIMAL:
    for indMouse, thisMouse in enumerate(np.unique(celldb.subject)):
        a = -1
        for indBinDV, thisQuantileDV in enumerate(quantilesDV):
            for indBinAP, thisQuantileAP in enumerate(quantilesAP):
                a = a+1
                quadrantTotalsByAnimal[indMouse, a] = np.sum(thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse))
                quadrantsSoundResponsiveByAnimal[indMouse, a] = np.sum(soundResponsive[thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse)])
                quadrantsSpeechResponsiveByAnimal[indMouse, a] = np.sum(speechResponsive[thisQuantileDV & thisQuantileAP & (celldb.subject == thisMouse)])
                quadrantsVotSelectiveByAnimal[indMouse, a] = np.sum(VOTselective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])
                quadrantsFtSelectiveByAnimal[indMouse, a] = np.sum(FTselective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])
                quadrantsMixedSelectiveByAnimal[indMouse, a] = np.sum(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])
                quadrantsSingleSelectiveByAnimal[indMouse, a] = np.sum(singleSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])


                #print(f'n Speech Selective {np.sum(singleSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)]) + np.sum(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])}')
                #print(f'n Vot Selective {np.sum(votSelective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])}')
                #print(f'n Ft Selective {np.sum(ftSelective[thisQuantileAP & thisQuantileDV & speechResponsive & (celldb.subject == thisMouse)])}')
                #print(f'n MixedSelective {np.sum(mixedSelective[thisQuantileDV & thisQuantileAP & speechResponsive & (celldb.subject == thisMouse)])}')




## -- run stats
## -- all cells
kstat, pValKruskalBestVOT_allcells = stats.kruskal(bestVotSIbyArea[0], bestVotSIbyArea[1], bestVotSIbyArea[2], bestVotSIbyArea[3], nan_policy = 'omit')
kstat, pValKruskalBestFT_allcells = stats.kruskal(bestFtSIbyArea[0], bestFtSIbyArea[1], bestFtSIbyArea[2], bestFtSIbyArea[3], nan_policy = 'omit')

# -- if group difference, test individual comparisons:
if pValKruskalBestVOT_allcells < 0.05:
    ustat, pValmannU_votAudPvsAudD_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[1])
    ustat, pValmannU_votAudPvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[2])
    ustat, pValmannU_votAudDvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[2])
    ustat, pValmannU_votAudPvsAudTea_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[3])
    ustat, pValmannU_votAudDvsAudTea_allcells = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[3])
    ustat, pValmannU_votAudVvsAudTea_allcells = stats.mannwhitneyu(bestVotSIbyArea[2], bestVotSIbyArea[3])

if pValKruskalBestFT_allcells < 0.05:
    ustat, pValmannU_ftAudPvsAudD_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[1])
    ustat, pValmannU_ftAudPvsAudV_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudDvsAudV_allcells = stats.mannwhitneyu(bestFtSIbyArea[1], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudPvsAudTea_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[3])
    ustat, pValmannU_ftAudDvsAudTea_allcells = stats.mannwhitneyu(bestFtSIbyArea[1], bestFtSIbyArea[3])
    ustat, pValmannU_ftAudVvsAudTea_allcells = stats.mannwhitneyu(bestFtSIbyArea[2], bestFtSIbyArea[3])

## -- just responsive cells
kstat, pValKruskalBestVOT = stats.kruskal(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2]], bestVotSIbyArea[3][speechResponsiveByArea[3]], nan_policy = 'omit')
kstat, pValKruskalBestFT = stats.kruskal(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2]], bestFtSIbyArea[3][speechResponsiveByArea[3]], nan_policy = 'omit')

## -- just selective cells
kstat, pValKruskalBestVOT_selective = stats.kruskal(bestVotSIbyArea[0][speechResponsiveByArea[0] & VOTselectivebyArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1] & VOTselectivebyArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2] & VOTselectivebyArea[2]], bestVotSIbyArea[3][speechResponsiveByArea[3] & VOTselectivebyArea[3]], nan_policy = 'omit')
kstat, pValKruskalBestFT_selective = stats.kruskal(bestFtSIbyArea[0][speechResponsiveByArea[0] & FTselectivebyArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1] & FTselectivebyArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2] & FTselectivebyArea[2]], bestFtSIbyArea[3][speechResponsiveByArea[3] & FTselectivebyArea[3]], nan_policy = 'omit')

# -- if group difference, test individual compariosns:
if pValKruskalBestVOT < 0.05:
    ustat, pValmannU_votAudPvsAudD = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_votAudPvsAudV = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_votAudDvsAudV = stats.mannwhitneyu(bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_votAudPvsTea = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[3][speechResponsiveByArea[3]])
    ustat, pValmannU_votAudDvsTea = stats.mannwhitneyu(bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[3][speechResponsiveByArea[3]])
    ustat, pValmannU_votAudVvsTea = stats.mannwhitneyu(bestVotSIbyArea[2][speechResponsiveByArea[2]], bestVotSIbyArea[3][speechResponsiveByArea[3]])

if pValKruskalBestFT < 0.05:
    ustat, pValmannU_ftAudPvsAudD = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_ftAudPvsAudV = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_ftAudDvsAudV = stats.mannwhitneyu(bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_ftAudPvsTea = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[3][speechResponsiveByArea[3]])
    ustat, pValmannU_ftAudDvsTea = stats.mannwhitneyu(bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[3][speechResponsiveByArea[3]])
    ustat, pValmannU_ftAudVvsTea = stats.mannwhitneyu(bestFtSIbyArea[2][speechResponsiveByArea[2]], bestFtSIbyArea[3][speechResponsiveByArea[3]])

## -- test proportion of speech responsive across cortical areas
oddsratio, pvalFracResponsive_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[1])],[np.sum(~speechResponsiveByArea[0]), np.sum(~speechResponsiveByArea[1])]]))
oddsratio, pvalFracResponsive_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[2])],[np.sum(~speechResponsiveByArea[0]), np.sum(~speechResponsiveByArea[2])]]))
oddsratio, pvalFracResponsive_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[2]), np.sum(speechResponsiveByArea[1])],[np.sum(~speechResponsiveByArea[2]), np.sum(~speechResponsiveByArea[1])]]))
oddsratio, pvalFracResponsive_AudPvsTea = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[3])],[np.sum(~speechResponsiveByArea[0]), np.sum(~speechResponsiveByArea[3])]]))
oddsratio, pvalFracResponsive_AudDvsTea = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[3]), np.sum(speechResponsiveByArea[1])],[np.sum(~speechResponsiveByArea[3]), np.sum(~speechResponsiveByArea[1])]]))
oddsratio, pvalFracResponsive_AudVvsTea = stats.fisher_exact(np.array([[np.sum(speechResponsiveByArea[3]), np.sum(speechResponsiveByArea[2])],[np.sum(~speechResponsiveByArea[3]), np.sum(~speechResponsiveByArea[2])]]))


## -- test selectivity distribution across cortical areas
oddsratio, pvalFracSelective_AudPvsAudD = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[1])]]))
oddsratio, pvalFracSelective_AudPvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[2])]]))
oddsratio, pvalFracSelective_AudDvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[2]), np.sum(notSelectivebyArea[1])]]))
oddsratio, pvalFracSelective_AudPvsTea = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[3]) + np.sum(mixedSelectivebyArea[3]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[3])]]))
oddsratio, pvalFracSelective_AudDvsTea = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1])), (np.sum(singleSelectivebyArea[3]) + np.sum(mixedSelectivebyArea[3]))],[np.sum(notSelectivebyArea[1]), np.sum(notSelectivebyArea[3])]]))
oddsratio, pvalFracSelective_AudVvsTea = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2])), (np.sum(singleSelectivebyArea[3]) + np.sum(mixedSelectivebyArea[3]))],[np.sum(notSelectivebyArea[2]), np.sum(notSelectivebyArea[3])]]))


## -- test distribution of feature selective cells b/w areas
oddsratio, pvalFracFtSelective_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(FTselectivebyArea[1][speechResponsiveByArea[1]])],[np.sum(~FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~FTselectivebyArea[1][speechResponsiveByArea[1]])]]))
oddsratio, pvalFracFtSelective_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(FTselectivebyArea[2][speechResponsiveByArea[2]])],[np.sum(~FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~FTselectivebyArea[2][speechResponsiveByArea[2]])]]))
oddsratio, pvalFracFtSelective_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(FTselectivebyArea[2][speechResponsiveByArea[2]])],[np.sum(~FTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(~FTselectivebyArea[2][speechResponsiveByArea[2]])]]))
oddsratio, pvalFracFtSelective_AudPvsTea = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(FTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~FTselectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracFtSelective_AudDvsTea = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(FTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(~FTselectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracFtSelective_AudVvsTea = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[2][speechResponsiveByArea[2]]), np.sum(FTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[2][speechResponsiveByArea[2]]), np.sum(~FTselectivebyArea[3][speechResponsiveByArea[3]])]]))


oddsratio, pvalFracVotSelective_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(VOTselectivebyArea[1][speechResponsiveByArea[1]])],[np.sum(~VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[1][speechResponsiveByArea[1]])]]))
oddsratio, pvalFracVotSelective_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(VOTselectivebyArea[2][speechResponsiveByArea[2]])],[np.sum(~VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[2][speechResponsiveByArea[2]])]]))
oddsratio, pvalFracVotSelective_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(VOTselectivebyArea[2][speechResponsiveByArea[2]])],[np.sum(~VOTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(~VOTselectivebyArea[2][speechResponsiveByArea[2]])]]))
oddsratio, pvalFracVotSelective_AudPvsTea = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(VOTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[0][speechResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracVotSelective_AudDvsTea = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(VOTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[1][speechResponsiveByArea[1]]), np.sum(~VOTselectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracVotSelective_AudVvsTea = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[2][speechResponsiveByArea[2]]), np.sum(VOTselectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[2][speechResponsiveByArea[2]]), np.sum(~VOTselectivebyArea[3][speechResponsiveByArea[3]])]]))


## -- test mixed selectivity
oddsratio, pvalFracMixed_AudPvsAudD = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(mixedSelectivebyArea[1][speechResponsiveByArea[1]])],[np.sum(singleSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(singleSelectivebyArea[1][speechResponsiveByArea[1]])]]))
oddsratio, pvalFracMixed_AudPvsAudV = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(mixedSelectivebyArea[2][speechResponsiveByArea[2]])],[np.sum(singleSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(singleSelectivebyArea[2][speechResponsiveByArea[2]])]]))
oddsratio, pvalFracMixed_AudDvsAudV = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[2][speechResponsiveByArea[2]]), np.sum(mixedSelectivebyArea[1][speechResponsiveByArea[1]])],[np.sum(singleSelectivebyArea[2][speechResponsiveByArea[2]]), np.sum(singleSelectivebyArea[1][speechResponsiveByArea[1]])]]))
oddsratio, pvalFracMixed_AudPvsTea = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(mixedSelectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(singleSelectivebyArea[0][speechResponsiveByArea[0]]), np.sum(singleSelectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracMixed_AudDvsTea = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[1][speechResponsiveByArea[1]]), np.sum(mixedSelectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(singleSelectivebyArea[1][speechResponsiveByArea[1]]), np.sum(singleSelectivebyArea[3][speechResponsiveByArea[3]])]]))
oddsratio, pvalFracMixed_AudVvsTea = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[2][speechResponsiveByArea[2]]), np.sum(mixedSelectivebyArea[3][speechResponsiveByArea[3]])],[np.sum(singleSelectivebyArea[2][speechResponsiveByArea[2]]), np.sum(singleSelectivebyArea[3][speechResponsiveByArea[3]])]]))



## -- test correlation b/w Speech feature selectivity and basic sound selectivity
kstat, pvalAmVot_AudP = stats.kruskal(bestVotSIbyArea[0][amSelectivebyArea[0] == 1],bestVotSIbyArea[0][amSelectivebyArea[0] == 0])
kstat, pvalAmVot_AudD = stats.kruskal(bestVotSIbyArea[1][amSelectivebyArea[1] == 1],bestVotSIbyArea[1][amSelectivebyArea[1] == 0])
kstat, pvalAmVot_AudV = stats.kruskal(bestVotSIbyArea[2][amSelectivebyArea[2] == 1],bestVotSIbyArea[2][amSelectivebyArea[2] == 0])
kstat, pvalAmVot_Tea = stats.kruskal(bestVotSIbyArea[3][amSelectivebyArea[3] == 1],bestVotSIbyArea[3][amSelectivebyArea[3] == 0])

kstat, pvalToneFt_AudP = stats.kruskal(bestFtSIbyArea[0][toneSelectivebyArea[0]==1], bestFtSIbyArea[0][toneSelectivebyArea[0]==0])
kstat, pvalToneFt_AudD = stats.kruskal(bestFtSIbyArea[1][toneSelectivebyArea[1]==1], bestFtSIbyArea[1][toneSelectivebyArea[1]==0])
kstat, pvalToneFt_AudV = stats.kruskal(bestFtSIbyArea[2][toneSelectivebyArea[2]==1], bestFtSIbyArea[2][toneSelectivebyArea[2]==0])
kstat, pvalToneFt_Tea = stats.kruskal(bestFtSIbyArea[3][toneSelectivebyArea[3]==1], bestFtSIbyArea[3][toneSelectivebyArea[3]==0])


## -- test if difference in proportion of selective cells among all cells (not just speech responsive)
oddsratio, pvalFracFtSelective_AudPvsAudD_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0]), np.sum(FTselectivebyArea[1])],[np.sum(~FTselectivebyArea[0]), np.sum(~FTselectivebyArea[1])]]))
oddsratio, pvalFracFtSelective_AudPvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0]), np.sum(FTselectivebyArea[2])],[np.sum(~FTselectivebyArea[0]), np.sum(~FTselectivebyArea[2])]]))
oddsratio, pvalFracFtSelective_AudDvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1]), np.sum(FTselectivebyArea[2])],[np.sum(~FTselectivebyArea[1]), np.sum(~FTselectivebyArea[2])]]))
oddsratio, pvalFracFtSelective_AudPvsTea_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0]), np.sum(FTselectivebyArea[3])],[np.sum(~FTselectivebyArea[0]), np.sum(~FTselectivebyArea[3])]]))
oddsratio, pvalFracFtSelective_AudDvsTea_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1]), np.sum(FTselectivebyArea[3])],[np.sum(~FTselectivebyArea[1]), np.sum(~FTselectivebyArea[3])]]))
oddsratio, pvalFracFtSelective_AudVvsTea_allcells = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[2]), np.sum(FTselectivebyArea[3])],[np.sum(~FTselectivebyArea[2]), np.sum(~FTselectivebyArea[3])]]))


oddsratio, pvalFracVotSelective_AudPvsAudD_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0]), np.sum(VOTselectivebyArea[1])],[np.sum(~VOTselectivebyArea[0]), np.sum(~VOTselectivebyArea[1])]]))
oddsratio, pvalFracVotSelective_AudPvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0]), np.sum(VOTselectivebyArea[2])],[np.sum(~VOTselectivebyArea[0]), np.sum(~VOTselectivebyArea[2])]]))
oddsratio, pvalFracVotSelective_AudDvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1]), np.sum(VOTselectivebyArea[2])],[np.sum(~VOTselectivebyArea[1]), np.sum(~VOTselectivebyArea[2])]]))
oddsratio, pvalFracVotSelective_AudPvsTea_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0]), np.sum(VOTselectivebyArea[3])],[np.sum(~VOTselectivebyArea[0]), np.sum(~VOTselectivebyArea[3])]]))
oddsratio, pvalFracVotSelective_AudDvsTea_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1]), np.sum(VOTselectivebyArea[3])],[np.sum(~VOTselectivebyArea[1]), np.sum(~VOTselectivebyArea[3])]]))
oddsratio, pvalFracVotSelective_AudVvsTea_allcells = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[2]), np.sum(VOTselectivebyArea[3])],[np.sum(~VOTselectivebyArea[2]), np.sum(~VOTselectivebyArea[3])]]))



## -- test if difference in proportion of selective cells among sound responsive cells (not just speech responsive)
oddsratio, pvalFracFtSelective_AudPvsAudD_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(FTselectivebyArea[1][soundResponsiveByArea[1]])],[np.sum(~FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~FTselectivebyArea[1][soundResponsiveByArea[1]])]]))
oddsratio, pvalFracFtSelective_AudPvsAudV_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(FTselectivebyArea[2][soundResponsiveByArea[2]])],[np.sum(~FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~FTselectivebyArea[2][soundResponsiveByArea[2]])]]))
oddsratio, pvalFracFtSelective_AudDvsAudV_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(FTselectivebyArea[2][soundResponsiveByArea[2]])],[np.sum(~FTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(~FTselectivebyArea[2][soundResponsiveByArea[2]])]]))
oddsratio, pvalFracFtSelective_AudPvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(FTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~FTselectivebyArea[3][soundResponsiveByArea[3]])]]))
oddsratio, pvalFracFtSelective_AudDvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(FTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(~FTselectivebyArea[3][soundResponsiveByArea[3]])]]))
oddsratio, pvalFracFtSelective_AudVvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(FTselectivebyArea[2][soundResponsiveByArea[2]]), np.sum(FTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~FTselectivebyArea[2][soundResponsiveByArea[2]]), np.sum(~FTselectivebyArea[3][soundResponsiveByArea[3]])]]))


oddsratio, pvalFracVotSelective_AudPvsAudD_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(VOTselectivebyArea[1][soundResponsiveByArea[1]])],[np.sum(~VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[1][soundResponsiveByArea[1]])]]))
oddsratio, pvalFracVotSelective_AudPvsAudV_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(VOTselectivebyArea[2][soundResponsiveByArea[2]])],[np.sum(~VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[2][soundResponsiveByArea[2]])]]))
oddsratio, pvalFracVotSelective_AudDvsAudV_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(VOTselectivebyArea[2][soundResponsiveByArea[2]])],[np.sum(~VOTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(~VOTselectivebyArea[2][soundResponsiveByArea[2]])]]))
oddsratio, pvalFracVotSelective_AudPvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(VOTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[0][soundResponsiveByArea[0]]), np.sum(~VOTselectivebyArea[3][soundResponsiveByArea[3]])]]))
oddsratio, pvalFracVotSelective_AudDvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(VOTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[1][soundResponsiveByArea[1]]), np.sum(~VOTselectivebyArea[3][soundResponsiveByArea[3]])]]))
oddsratio, pvalFracVotSelective_AudVvsTea_soundResponsive = stats.fisher_exact(np.array([[np.sum(VOTselectivebyArea[2][soundResponsiveByArea[2]]), np.sum(VOTselectivebyArea[3][soundResponsiveByArea[3]])],[np.sum(~VOTselectivebyArea[2][soundResponsiveByArea[2]]), np.sum(~VOTselectivebyArea[3][soundResponsiveByArea[3]])]]))





if STATSUMMARY:
    print('--Stats Summary--')
    print(f'Bonferroni corrected alpha = {np.round(0.05/6, 3)}')
    print(f'Min FR threshold for best SI = {frMinBestSI}sp/s')
    print(f'Group differences SI-VOT: p = {np.round(pValKruskalBestVOT_selective,3)}')
    print(f'Group differences SI-FT: p = {np.round(pValKruskalBestFT_selective,3)}')
    if pValKruskalBestVOT_selective < 0.05:
        print(f'SI-VOT AudP vs AudD: p = {np.round(pValmannU_votAudPvsAudD,3)}')
        print(f'SI-VOT AudP vs AudV: p = {np.round(pValmannU_votAudPvsAudV,3)}')
        print(f'SI-VOT AudD vs AudV: p = {np.round(pValmannU_votAudDvsAudV,3)}')
        print(f'SI-VOT AudP vs. Tea: p = {np.round(pValmannU_votAudPvsTea,3)}')
        print(f'SI-VOT AudD vs. Tea: p = {np.round(pValmannU_votAudDvsTea,3)}')
        print(f'SI-VOT AudV vs. Tea: p = {np.round(pValmannU_votAudVvsTea,3)}')
    if pValKruskalBestFT_selective < 0.05:
        print(f'SI-FT AudP vs AudD: p = {np.round(pValmannU_ftAudPvsAudD,3)}')
        print(f'SI-FT AudP vs AudV: p = {np.round(pValmannU_ftAudPvsAudV,3)}')
        print(f'SI-FT AudD vs AudV: p = {np.round(pValmannU_ftAudDvsAudV,3)}')
        print(f'SI-FT AudP vs. Tea: p = {np.round(pValmannU_ftAudPvsTea,3)}')
        print(f'SI-FT AudD vs. Tea: p = {np.round(pValmannU_ftAudDvsTea,3)}')
        print(f'SI-FT AudV vs. Tea: p = {np.round(pValmannU_ftAudVvsTea,3)}')
    print(f'Frac VOT selective AudP vs AudD p = {np.round(pvalFracVotSelective_AudPvsAudD,3)}')
    print(f'Frac VOT selective AudP vs AudV p = {np.round(pvalFracVotSelective_AudPvsAudV,3)}')
    print(f'Frac VOT selective AudD vs AudV p = {np.round(pvalFracVotSelective_AudDvsAudV,3)}')
    print(f'Frac VOT selective AudP vs Tea p = {np.round(pvalFracVotSelective_AudPvsTea,3)}')
    print(f'Frac VOT selective AudD vs Tea p = {np.round(pvalFracVotSelective_AudDvsTea,3)}')
    print(f'Frac VOT selective AudV vs Tea p = {np.round(pvalFracVotSelective_AudVvsTea,3)}')
    print(f'Frac FT selective AudP vs AudD p = {np.round(pvalFracFtSelective_AudPvsAudD,3)}')
    print(f'Frac FT selective AudP vs AudV p = {np.round(pvalFracFtSelective_AudPvsAudV,3)}')
    print(f'Frac FT selective AudD vs AudV p = {np.round(pvalFracFtSelective_AudDvsAudV,3)}')
    print(f'Frac FT selective AudP vs Tea p = {np.round(pvalFracFtSelective_AudPvsTea,3)}')
    print(f'Frac FT selective AudD vs Tea p = {np.round(pvalFracFtSelective_AudDvsTea,3)}')
    print(f'Frac FT selective AudV vs Tea p = {np.round(pvalFracFtSelective_AudVvsTea,3)}')
    print(f'Frac mixed selective AudP vs AudD:  p = {np.round(pvalFracMixed_AudPvsAudD,3)}')
    print(f'Frac mixed selective AudP vs AudV:  p = {np.round(pvalFracMixed_AudPvsAudV,3)}')
    print(f'Frac mixed selective AudD vs AudV:  p = {np.round(pvalFracMixed_AudDvsAudV,3)}')
    print(f'Frac mixed selective AudP vs Tea:  p = {np.round(pvalFracMixed_AudPvsTea,3)}')
    print(f'Frac mixed selective AudD vs Tea:  p = {np.round(pvalFracMixed_AudDvsTea,3)}')
    print(f'Frac mixed selective AudV vs Tea:  p = {np.round(pvalFracMixed_AudVvsTea,3)}')
    print(f'AudP n: {len(speechResponsiveByArea[0])}, n speechResponsive: {np.sum(speechResponsiveByArea[0])}, n selective: VOT = {np.sum(VOTselectivebyArea[0][speechResponsiveByArea[0]])} ({np.round(np.mean(VOTselectivebyArea[0][speechResponsiveByArea[0]])*100,1)}%), FT = {np.sum(FTselectivebyArea[0][speechResponsiveByArea[0]])} ({np.round(np.mean(FTselectivebyArea[0][speechResponsiveByArea[0]])*100,1)}%), Mixed = {np.sum(mixedSelectivebyArea[0][speechResponsiveByArea[0]])} ({np.round(np.mean(mixedSelectivebyArea[0][speechResponsiveByArea[0]])*100,1)}%)')
    print(f'AudD n: {len(speechResponsiveByArea[1])}, n speechResponsive: {np.sum(speechResponsiveByArea[1][speechResponsiveByArea[1]])}, n selective: VOT = {np.sum(VOTselectivebyArea[1][speechResponsiveByArea[1]])} ({np.round(np.mean(VOTselectivebyArea[1][speechResponsiveByArea[1]])*100,1)}%), FT = {np.sum(FTselectivebyArea[1][speechResponsiveByArea[1]])} ({np.round(np.mean(FTselectivebyArea[1][speechResponsiveByArea[1]])*100,1)}%), Mixed = {np.sum(mixedSelectivebyArea[1][speechResponsiveByArea[1]])} ({np.round(np.mean(mixedSelectivebyArea[1][speechResponsiveByArea[1]])*100,1)}%)')
    print(f'AudV n: {len(speechResponsiveByArea[2])}, n speechResponsive: {np.sum(speechResponsiveByArea[2])}, n selective: VOT = {np.sum(VOTselectivebyArea[2][speechResponsiveByArea[2]])} ({np.round(np.mean(VOTselectivebyArea[2][speechResponsiveByArea[2]])*100,1)}%), FT = {np.sum(FTselectivebyArea[2][speechResponsiveByArea[2]])} ({np.round(np.mean(FTselectivebyArea[2][speechResponsiveByArea[2]])*100,1)}%), Mixed = {np.sum(mixedSelectivebyArea[2][speechResponsiveByArea[2]])} ({np.round(np.mean(mixedSelectivebyArea[2][speechResponsiveByArea[2]])*100,1)}%)')
    print(f'Tea n: {len(speechResponsiveByArea[3])}, n speechResponsive: {np.sum(speechResponsiveByArea[3])}, n selective: VOT = {np.sum(VOTselectivebyArea[3][speechResponsiveByArea[3]])} ({np.round(np.mean(VOTselectivebyArea[3][speechResponsiveByArea[3]])*100,1)}%), FT = {np.sum(FTselectivebyArea[3][speechResponsiveByArea[3]])} ({np.round(np.mean(FTselectivebyArea[3][speechResponsiveByArea[3]])*100,1)}%), Mixed = {np.sum(mixedSelectivebyArea[3][speechResponsiveByArea[3]])} ({np.round(np.mean(mixedSelectivebyArea[3][speechResponsiveByArea[3]])*100,1)}%)')




np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityIndexFt = bestSelectivityIndexFt, bestSelectivityIndexVot = bestSelectivityIndexVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, excludeSpeech = excludeSpeech, pValKruskalBestFT = pValKruskalBestFT, pValKruskalBestVOT = pValKruskalBestVOT, speechResponsive = speechResponsive, amResponsive = amResponsive, toneResponsive = toneResponsive, soundResponsive = soundResponsive, amSelective = amSelective, toneSelective = toneSelective, maxFiringRateSpeechEvoked = maxFiringRateSpeechEvoked, isAudArea = isAudArea, y_coord = y_coords, z_coord = z_coords, x_coord = x_coords, z_coords_jittered = z_coords_jittered, x_coords_jittered = x_coords_jittered, subject = celldb.subject, date = celldb.date, cluster = celldb.cluster, pvalPermutationtestFt = pvalPermutationtestFt, pvalPermutationtestVot = pvalPermutationtestVot, shuffledVotBest = shuffledVotBest, shuffledFtBest = shuffledFtBest, whichFT = whichFT, whichVOT = whichVOT, whichFT_labels = whichFT_labels, whichVOT_labels = whichVOT_labels, isCortical = isCortical, quantilesDV = quantilesDV, quantilesAP = quantilesAP, quadrantBoundsDV = quadrantBoundsDV, quadrantBoundsAP = quadrantBoundsAP, quadrantBoundsDV_AApix = quadrantBoundsDV_AApix, quadrantBoundsAP_AApix = quadrantBoundsAP_AApix, quadrantLabels = quadrantLabels, quadrantTotals = quadrantTotals, quadrantsVOT = quadrantsVOT, quadrantsFT = quadrantsFT, quadrantsVotSelective = quadrantsVotSelective, quadrantsFtSelective = quadrantsFtSelective, quadrantsMixedSelective = quadrantsMixedSelective, quadrantsSingleSelective = quadrantsSingleSelective, quadrantsSpeechResponsive = quadrantsSpeechResponsive, quadrantsSoundResponsive = quadrantsSoundResponsive, quadrantTotalsByAnimal = quadrantTotalsByAnimal, quadrantsSoundResponsiveByAnimal = quadrantsSoundResponsiveByAnimal, quadrantsSpeechResponsiveByAnimal = quadrantsSpeechResponsiveByAnimal, quadrantsVotSelectiveByAnimal = quadrantsVotSelectiveByAnimal, quadrantsFtSelectiveByAnimal = quadrantsFtSelectiveByAnimal, quadrantsSingleSelectiveByAnimal = quadrantsSingleSelectiveByAnimal, quadrantsMixedSelectiveByAnimal = quadrantsMixedSelectiveByAnimal, layersDeep = layersDeep, layer4 = layer4, layersSuperficial = layersSuperficial )
print('saved to ' f'{figDataFullPath}')
