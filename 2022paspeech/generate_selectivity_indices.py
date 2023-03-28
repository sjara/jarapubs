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
STATSUMMARY = 0
AREASBYDEPTH = 1

if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)

figDataFullPath = os.path.join(figDataDir,figDataFile)
#scriptFullPath = os.path.realpath(__file__)

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)
isAudArea = celldb.recordingAreaName.str.contains('auditory area')


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

randSI = np.zeros(10000)
for indPerm, thisPerm in enumerate(randSI):
    randFRs = np.random.randint(1,60,2)
    randSI[indPerm] = (np.max(randFRs) - np.min(randFRs))/(np.max(randFRs) + np.min(randFRs))


alpha = 0.05/12
speechResponsive = np.array(celldb.speechMinPvalOnset < alpha)

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

## -- set selectivityCriterion for mixed selectivity
selectivityCriterion = 0.4
VOTselective = bestSelectivityIndexVot > selectivityCriterion
FTselective = bestSelectivityIndexFt > selectivityCriterion
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
    '''
    audP_coords = (y_coords > 106) & (y_coords <= 124) #& isAudArea #& (celldb.z_coord > 185)
    audD_coords = (y_coords <= 106) #& isAudArea #& (celldb.z_coord > 185)
    audV_coords = (y_coords > 124) #& isAudArea #& (celldb.z_coord > 185)
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

    audAreas_byCoords = np.array([audP_coords, audD_coords, audV_coords])


    for indArea, thisArea in enumerate(audCtxAreas):
        bestFtSIbyArea.append(bestSelectivityIndexFt[audAreas_byCoords[indArea]])
        bestVotSIbyArea.append(bestSelectivityIndexVot[audAreas_byCoords[indArea]])
        speechResponsiveByArea.append(speechResponsive[audAreas_byCoords[indArea]])
        amSelectivebyArea.append(amSelective[audAreas_byCoords[indArea]])
        toneSelectivebyArea.append(toneSelective[audAreas_byCoords[indArea]])
        excludeCellsbyArea.append(excludeCells[audAreas_byCoords[indArea]])
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
    mixedSelectivebyArea[indArea] = mixedSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    singleSelectivebyArea[indArea] = singleSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    notSelectivebyArea[indArea] = notSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    y_coordsbyArea[indArea] = y_coordsbyArea[indArea][~excludeCellsbyArea[indArea]]
    z_coordsbyArea[indArea] = z_coordsbyArea[indArea][~excludeCellsbyArea[indArea]]




## -- run stats
## -- all cells
kstat, pValKruskalBestVOT_allcells = stats.kruskal(bestVotSIbyArea[0], bestVotSIbyArea[1], bestVotSIbyArea[2], nan_policy = 'omit')
kstat, pValKruskalBestFT_allcells = stats.kruskal(bestFtSIbyArea[0], bestFtSIbyArea[1], bestFtSIbyArea[2], nan_policy = 'omit')

# -- if group difference, test individual comparisons:
if pValKruskalBestVOT_allcells < 0.05:
    ustat, pValmannU_votAudPvsAudD_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[1])
    ustat, pValmannU_votAudPvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[2])
    ustat, pValmannU_votAudDvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[2])

if pValKruskalBestFT_allcells < 0.05:
    ustat, pValmannU_ftAudPvsAudD_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[1])
    ustat, pValmannU_ftAudPvsAudV_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudDvsAudV_allcells = stats.mannwhitneyu(bestFtSIbyArea[1], bestFtSIbyArea[2])

## -- just responsive cells
kstat, pValKruskalBestVOT = stats.kruskal(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2]], nan_policy = 'omit')
kstat, pValKruskalBestFT = stats.kruskal(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2]], nan_policy = 'omit')

# -- if group difference, test individual compariosns:
if pValKruskalBestVOT < 0.05:
    ustat, pValmannU_votAudPvsAudD = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_votAudPvsAudV = stats.mannwhitneyu(bestVotSIbyArea[0][speechResponsiveByArea[0]], bestVotSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_votAudDvsAudV = stats.mannwhitneyu(bestVotSIbyArea[1][speechResponsiveByArea[1]], bestVotSIbyArea[2][speechResponsiveByArea[2]])

if pValKruskalBestFT < 0.05:
    ustat, pValmannU_ftAudPvsAudD = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[1][speechResponsiveByArea[1]])
    ustat, pValmannU_ftAudPvsAudV = stats.mannwhitneyu(bestFtSIbyArea[0][speechResponsiveByArea[0]], bestFtSIbyArea[2][speechResponsiveByArea[2]])
    ustat, pValmannU_ftAudDvsAudV = stats.mannwhitneyu(bestFtSIbyArea[1][speechResponsiveByArea[1]], bestFtSIbyArea[2][speechResponsiveByArea[2]])

## -- test selectivity distribution across cortical areas
oddsratio, pvalFracSelective_AudPvsAudD = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[1])]]))
oddsratio, pvalFracSelective_AudPvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[2])]]))
oddsratio, pvalFracSelective_AudDvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[2]), np.sum(notSelectivebyArea[1])]]))

## -- test mixed selectivity
oddsratio, pvalFracMixed_AudPvsAudD = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[1])]]))
oddsratio, pvalFracMixed_AudPvsAudV = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[2])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[2])]]))
oddsratio, pvalFracMixed_AudDvsAudV = stats.fisher_exact(2*np.array([[np.sum(mixedSelectivebyArea[2]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[2]), np.sum(singleSelectivebyArea[1])]]))

## -- test correlation b/w Speech feature selectivity and basic sound selectivity
kstat, pvalAmVot_AudP = stats.kruskal(bestVotSIbyArea[0][amSelectivebyArea[0] == 1],bestVotSIbyArea[0][amSelectivebyArea[0] == 0])
kstat, pvalAmVot_AudD = stats.kruskal(bestVotSIbyArea[1][amSelectivebyArea[1] == 1],bestVotSIbyArea[1][amSelectivebyArea[1] == 0])
kstat, pvalAmVot_AudV = stats.kruskal(bestVotSIbyArea[2][amSelectivebyArea[2] == 1],bestVotSIbyArea[2][amSelectivebyArea[2] == 0])

kstat, pvalToneFt_AudP = stats.kruskal(bestFtSIbyArea[0][toneSelectivebyArea[0]==1], bestFtSIbyArea[0][toneSelectivebyArea[0]==0])
kstat, pvalToneFt_AudD = stats.kruskal(bestFtSIbyArea[1][toneSelectivebyArea[1]==1], bestFtSIbyArea[1][toneSelectivebyArea[1]==0])
kstat, pvalToneFt_AudV = stats.kruskal(bestFtSIbyArea[2][toneSelectivebyArea[2]==1], bestFtSIbyArea[2][toneSelectivebyArea[2]==0])



if STATSUMMARY:
    print('--Stats Summary--')
    print(f'Group differences SI-FT: p = {np.round(pValKruskalBestFT,3)}')
    print(f'Group differences SI-VOT: p = {np.round(pValKruskalBestVOT,3)}')
    if pValKruskalBestVOT < 0.05:
        print(f'SI-VOT AudP vs AudD: p = {np.round(pValmannU_votAudPvsAudD,3)}')
        print(f'SI-VOT AudP vs AudV: p = {np.round(pValmannU_votAudPvsAudV,3)}')
        print(f'SI-VOT AudD vs AudV: p = {np.round(pValmannU_votAudDvsAudV,3)}')
    if pValKruskalBestFT < 0.05:
        print(f'SI-FT AudP vs AudD: p = {np.round(pValmannU_ftAudPvsAudD,3)}')
        print(f'SI-FT AudP vs AudV: p = {np.round(pValmannU_ftAudPvsAudV,3)}')
        print(f'SI-FT AudD vs AudV: p = {np.round(pValmannU_ftAudDvsAudV,3)}')
    print(f'exclusion criterion = firing rate < {exclusionCriterion} sp/s')
    print(f'n Excluded {audCtxAreas[0]}: {np.sum(excludeCellsbyArea[0])}')
    print(f'n Excluded {audCtxAreas[1]}: {np.sum(excludeCellsbyArea[1])}')
    print(f'n Excluded {audCtxAreas[2]}: {np.sum(excludeCellsbyArea[2])}')
    print(f'Selectivity criterion = {selectivityCriterion}')
    #print(f'Frac selective AudP vs AudD: p = {pvalFracSelective_AudPvsAudD}')
    #print(f'Frac selective AudP vs AudV: p = {pvalFracSelective_AudPvsAudV}')
    #print(f'Frac selective AudD vs AudV: p = {pvalFracSelective_AudDvsAudV}')
    print(f'Frac mixed selective AudP vs AudD:  p = {pvalFracMixed_AudPvsAudD}')
    print(f'Frac mixed selective AudP vs AudV:  p = {pvalFracMixed_AudPvsAudV}')
    print(f'Frac mixed selective AudD vs AudV:  p = {pvalFracMixed_AudDvsAudV}')
    print(f'Frac Selective AudP = {(np.sum(mixedSelectivebyArea[0]) + np.sum(singleSelectivebyArea[0]))/len(singleSelectivebyArea[0])}')
    print(f'Frac Selective AudD = {(np.sum(mixedSelectivebyArea[1]) + np.sum(singleSelectivebyArea[1]))/len(singleSelectivebyArea[1])}')
    print(f'Frac Selective AudV = {(np.sum(mixedSelectivebyArea[2]) + np.sum(singleSelectivebyArea[2]))/len(singleSelectivebyArea[2])}')
    print(f'AudP n speechResponsive {np.sum(speechResponsiveByArea[0])}, n selective: VOT = {np.sum(bestVotSIbyArea[0][speechResponsiveByArea[0]]>selectivityCriterion)} ({np.round((np.sum(bestVotSIbyArea[0][speechResponsiveByArea[0]]>selectivityCriterion)/np.sum(speechResponsiveByArea[0]))*100,1)}%), FT = {np.sum(bestFtSIbyArea[0][speechResponsiveByArea[0]]>selectivityCriterion)} ({np.round((np.sum(bestFtSIbyArea[0][speechResponsiveByArea[0]]>selectivityCriterion)/np.sum(speechResponsiveByArea[0]))*100,1)}%)')
    print(f'AudD n speechResponsive {np.sum(speechResponsiveByArea[1])}, n selective: VOT = {np.sum(bestVotSIbyArea[1][speechResponsiveByArea[1]]>selectivityCriterion)} ({np.round((np.sum(bestVotSIbyArea[1][speechResponsiveByArea[1]]>selectivityCriterion)/np.sum(speechResponsiveByArea[1]))*100,1)}%), FT = {np.sum(bestFtSIbyArea[1][speechResponsiveByArea[1]]>selectivityCriterion)} ({np.round((np.sum(bestFtSIbyArea[1][speechResponsiveByArea[1]]>selectivityCriterion)/np.sum(speechResponsiveByArea[1]))*100,1)}%)')
    print(f'AudV n speechResponsive {np.sum(speechResponsiveByArea[2])}, n selective: VOT = {np.sum(bestVotSIbyArea[2][speechResponsiveByArea[2]]>selectivityCriterion)} ({np.round((np.sum(bestVotSIbyArea[2][speechResponsiveByArea[2]]>selectivityCriterion)/np.sum(speechResponsiveByArea[2]))*100,1)}%), FT = {np.sum(bestFtSIbyArea[2][speechResponsiveByArea[2]]>selectivityCriterion)} ({np.round((np.sum(bestFtSIbyArea[2][speechResponsiveByArea[2]]>selectivityCriterion)/np.sum(speechResponsiveByArea[2]))*100,1)}%)')

np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityIndexFt = bestSelectivityIndexFt, bestSelectivityIndexVot = bestSelectivityIndexVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, excludeCells = excludeCells, pValKruskalBestFT = pValKruskalBestFT, pValKruskalBestVOT = pValKruskalBestVOT, speechResponsive = speechResponsive, amResponsive = amResponsive, toneResponsive = toneResponsive, soundResponsive = soundResponsive, amSelective = amSelective, toneSelective = toneSelective, maxFiringRateSpeechEvoked = maxFiringRateSpeechEvoked, isAudArea = isAudArea, y_coord = y_coords, z_coord = z_coords, x_coord = x_coords, z_coords_jittered = z_coords_jittered, x_coords_jittered = x_coords_jittered, randSI = randSI)
print('saved to ' f'{figDataFullPath}')
