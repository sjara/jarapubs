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

FIGNAME = 'selectivityIndices'
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
STATSUMMARY = 1
AREASBYDEPTH = 0
shuffledDataFile = 'data_shuffledSIs.npz'
#figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)


if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)

figDataFullPath = os.path.join(figDataDir,figDataFile)
shuffledDataFullPath = os.path.join(figDataDir,shuffledDataFile)
shuffledData = np.load(shuffledDataFullPath, allow_pickle=True)
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
ctxAreas = ['Anterolateral visual area', 'Dorsal auditory area', 'Ectorhinal area', 'Endopiriform nucleus', 'Lateral visual area', 'Laterointermediate area', 'Perirhinal area', 'Posterior auditory area', 'Primary auditory area', 'Supplemental somatosensory area', 'Temporal association areas', 'Ventral auditory area']
isCortical = np.zeros(len(isAudArea), dtype = bool)
for indArea, thisArea in enumerate(ctxAreas):
    isCortical[recordingAreaName == thisArea] = True


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

speechAlpha = 0.05/12
speechResponsive = np.array(celldb.speechMinPvalOnset < speechAlpha)
amAlpha = 0.05/11 #11 AM frequencies presented
toneAlpha = 0.05/16 #16 pure tones presented
toneResponsive = np.array(celldb.toneMinPval < toneAlpha)
amResponsive = np.array((celldb.amMinPvalOnset < amAlpha) |
                    (celldb.amMinPvalSustain < amAlpha))
soundResponsive = np.array(toneResponsive|amResponsive|speechResponsive)


## -- exclude cells with low firing rates
exclusionCriterion = 3  # sp/s
maxFiringRateSpeechEvoked = np.max([celldb.maxFiringRate_FT_VOTmin, celldb.maxFiringRate_FT_VOTmax, celldb.maxFiringRate_VOT_FTmin, celldb.maxFiringRate_VOT_FTmax], 0)
exclude_evoked = maxFiringRateSpeechEvoked < exclusionCriterion
exclude_baseline = celldb.speechFiringRateBaseline < exclusionCriterion
excludeCells = exclude_evoked & exclude_baseline

## -- select best index b/w conditions
ftCombos = np.array([selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax])
votCombos = np.array([selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax])
bestSelectivityIndexFt = ftCombos.max(0)
bestSelectivityIndexVot = votCombos.max(0)

amOnsetSelective = np.array(celldb.amSelectivityPvalOnset < 0.05)
amSustainSelective = np.array(celldb.amSelectivityPvalSustain < 0.05)
amSelective = amOnsetSelective| amSustainSelective
toneSelective = np.array(celldb.toneSelectivityPval < 0.05)
excludeAM = celldb.amFiringRateMaxOnset < exclusionCriterion
excludeTone = celldb.toneFiringRateBest < exclusionCriterion
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

pvalPermutationtestFt = np.ones(nCells)
pvalPermutationtestVot = np.ones(nCells)
for indCell, thisCell in enumerate(votSI_FTmax):
    if votSI_FTmax[indCell] > votSI_FTmin[indCell]:
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmax[indCell,:] >=votSI_FTmax[indCell])
    elif votSI_FTmax[indCell] < votSI_FTmin[indCell]:
        pvalPermutationtestVot[indCell] = np.mean(shuffledSIVot_FTmin[indCell,:] >=votSI_FTmin[indCell])

    if ftSI_VOTmax[indCell] > ftSI_VOTmin[indCell]:
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmax[indCell,:] >=ftSI_VOTmax[indCell])
    elif ftSI_VOTmax[indCell] < ftSI_VOTmin[indCell]:
        pvalPermutationtestFt[indCell] = np.mean(shuffledSIFt_VOTmin[indCell,:] >= ftSI_VOTmin[indCell])


'''
pvalPermutationtestFt = np.ones(len(bestSelectivityIndexFt))
pvalPermutationtestVot = np.ones(len(bestSelectivityIndexFt))

for indCell, thisCell in enumerate(bestSelectivityIndexFt):
    pvalPermutationtestFt[indCell] = np.mean(shuffledSiEachCell[indCell,:] >= bestSelectivityIndexFt[indCell])
    pvalPermutationtestVot[indCell] = np.mean(shuffledSiEachCell[indCell,:] >= bestSelectivityIndexVot[indCell])
'''

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
speechResponsiveByArea = []
VOTselectivebyArea = []
FTselectivebyArea = []
amSelectivebyArea = []
toneSelectivebyArea = []
excludeCellsbyArea = []
mixedSelectivebyArea = []
singleSelectivebyArea = []
notSelectivebyArea = []
y_coordsbyArea = []
z_coordsbyArea = []

if AREASBYDEPTH == 0:
    for indArea, thisArea in enumerate(audCtxAreas):
        bestFtSIbyArea.append(bestSelectivityIndexFt[recordingAreaName == thisArea])
        bestVotSIbyArea.append(bestSelectivityIndexVot[recordingAreaName == thisArea])
        speechResponsiveByArea.append(speechResponsive[recordingAreaName == thisArea])
        amSelectivebyArea.append(amSelective[recordingAreaName == thisArea])
        toneSelectivebyArea.append(toneSelective[recordingAreaName == thisArea])
        excludeCellsbyArea.append(excludeCells[recordingAreaName == thisArea])
        VOTselectivebyArea.append(VOTselective[recordingAreaName == thisArea])
        FTselectivebyArea.append(FTselective[recordingAreaName == thisArea])
        mixedSelectivebyArea.append(mixedSelective[recordingAreaName == thisArea])
        singleSelectivebyArea.append(singleSelective[recordingAreaName == thisArea])
        notSelectivebyArea.append(notSelective[recordingAreaName == thisArea])
        y_coordsbyArea.append(y_coords[recordingAreaName == thisArea])
        z_coordsbyArea.append(z_coords[recordingAreaName == thisArea])

if AREASBYDEPTH: ## used for test where we defined ACtx area by depth
    print('using depth to define auditory areas')
    isAudArea = celldb.recordingAreaName.str.contains('auditory area')
    # Coords as seen as splits on hist: audD < 105 < audP < 133 < audV
    # Coords as even thirds: audD < 110 < audP < 124 < audV

    audP_coords = (y_coords > 106) & (y_coords <= 124) & isCortical#& isAudArea #& (celldb.z_coord > 185)
    audD_coords = (y_coords <= 106) & isCortical #& isAudArea #& (celldb.z_coord > 185)
    audV_coords = (y_coords > 124) & (y_coords <=145) & isCortical#& isAudArea #& (celldb.z_coord > 185)
    tea_coords = (y_coords > 145) & isCortical
    '''
    nBins = 4
    binSize = (np.max(y_coords[speechResponsive & ~excludeCells]) - np.min(y_coords[speechResponsive & ~excludeCells]))/nBins
    bins = np.arange(np.min(y_coords[speechResponsive & ~excludeCells]), np.max(y_coords[speechResponsive & ~excludeCells]), binSize)
    audD_coords = (y_coords >= bins[0]) & (y_coords < bins[1])
    audP_coords = (y_coords >= bins[1]) & (y_coords < bins[2])
    if nBins == 3:
        audV_coords = y_coords >= bins[2]
    elif nBins == 4:
        audV_coords = (y_coords >= bins[2]) & (y_coords < bins[3])
    '''
    audAreas_byCoords = np.array([audP_coords, audD_coords, audV_coords, tea_coords])


    for indArea, thisArea in enumerate(audAreas_byCoords):
        bestFtSIbyArea.append(bestSelectivityIndexFt[audAreas_byCoords[indArea]])
        bestVotSIbyArea.append(bestSelectivityIndexVot[audAreas_byCoords[indArea]])
        speechResponsiveByArea.append(speechResponsive[audAreas_byCoords[indArea]])
        amSelectivebyArea.append(amSelective[audAreas_byCoords[indArea]])
        toneSelectivebyArea.append(toneSelective[audAreas_byCoords[indArea]])
        excludeCellsbyArea.append(excludeCells[audAreas_byCoords[indArea]])
        VOTselectivebyArea.append(VOTselective[audAreas_byCoords[indArea]])
        FTselectivebyArea.append(FTselective[audAreas_byCoords[indArea]])
        mixedSelectivebyArea.append(mixedSelective[audAreas_byCoords[indArea]])
        singleSelectivebyArea.append(singleSelective[audAreas_byCoords[indArea]])
        notSelectivebyArea.append(notSelective[audAreas_byCoords[indArea]])
        y_coordsbyArea.append(y_coords[audAreas_byCoords[indArea]])
        z_coordsbyArea.append(z_coords[audAreas_byCoords[indArea]])


## -- exclude low spike count cells
for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea[indArea] = bestFtSIbyArea[indArea][~excludeCellsbyArea[indArea]]
    bestVotSIbyArea[indArea] = bestVotSIbyArea[indArea][~excludeCellsbyArea[indArea]]
    speechResponsiveByArea[indArea] = speechResponsiveByArea[indArea][~excludeCellsbyArea[indArea]]
    amSelectivebyArea[indArea] = amSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    toneSelectivebyArea[indArea] = toneSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    VOTselectivebyArea[indArea] = VOTselectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    FTselectivebyArea[indArea] = FTselectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    mixedSelectivebyArea[indArea] = mixedSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    singleSelectivebyArea[indArea] = singleSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    notSelectivebyArea[indArea] = notSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    y_coordsbyArea[indArea] = y_coordsbyArea[indArea][~excludeCellsbyArea[indArea]]
    z_coordsbyArea[indArea] = z_coordsbyArea[indArea][~excludeCellsbyArea[indArea]]



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



if STATSUMMARY:
    print('--Stats Summary--')
    print(f'Group differences SI-VOT: p = {np.round(pValKruskalBestVOT,3)}')
    print(f'Group differences SI-FT: p = {np.round(pValKruskalBestFT,3)}')
    print(f'Bonferroni corrected alpha = {np.round(0.05/6, 3)}')
    if pValKruskalBestVOT < 0.05:
        print(f'SI-VOT AudP vs AudD: p = {np.round(pValmannU_votAudPvsAudD,3)}')
        print(f'SI-VOT AudP vs AudV: p = {np.round(pValmannU_votAudPvsAudV,3)}')
        print(f'SI-VOT AudD vs AudV: p = {np.round(pValmannU_votAudDvsAudV,3)}')
        print(f'SI-VOT AudP vs. Tea: p = {np.round(pValmannU_votAudPvsTea,3)}')
        print(f'SI-VOT AudD vs. Tea: p = {np.round(pValmannU_votAudDvsTea,3)}')
        print(f'SI-VOT AudV vs. Tea: p = {np.round(pValmannU_votAudVvsTea,3)}')
    if pValKruskalBestFT < 0.05:
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
    print(f'exclusion criterion = firing rate < {exclusionCriterion} sp/s')
    print(f'n Excluded {audCtxAreas[0]}: {np.sum(excludeCellsbyArea[0])}')
    print(f'n Excluded {audCtxAreas[1]}: {np.sum(excludeCellsbyArea[1])}')
    print(f'n Excluded {audCtxAreas[2]}: {np.sum(excludeCellsbyArea[2])}')
    print(f'n Excluded {audCtxAreas[3]}: {np.sum(excludeCellsbyArea[3])}')
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

'''
np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityIndexFt = bestSelectivityIndexFt, bestSelectivityIndexVot = bestSelectivityIndexVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, excludeCells = excludeCells, pValKruskalBestFT = pValKruskalBestFT, pValKruskalBestVOT = pValKruskalBestVOT, speechResponsive = speechResponsive, amResponsive = amResponsive, toneResponsive = toneResponsive, soundResponsive = soundResponsive, amSelective = amSelective, toneSelective = toneSelective, maxFiringRateSpeechEvoked = maxFiringRateSpeechEvoked, isAudArea = isAudArea, y_coord = y_coords, z_coord = z_coords, x_coord = x_coords, z_coords_jittered = z_coords_jittered, x_coords_jittered = x_coords_jittered, subject = celldb.subject, date = celldb.date, cluster = celldb.cluster, pvalPermutationtestFt = pvalPermutationtestFt, pvalPermutationtestVot = pvalPermutationtestVot)
print('saved to ' f'{figDataFullPath}')
'''
