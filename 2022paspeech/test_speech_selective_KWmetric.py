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

FIGNAME = 'figure_neuropix_ac'
figDataFile = 'data_selectivity_indices.npz'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
STATSUMMARY = 1
AREASBYDEPTH = 1

if not os.path.exists(figDataDir):
    os.mkdir(figDataDir)

figDataFullPath = os.path.join(figDataDir,figDataFile)
scriptFullPath = os.path.realpath(__file__)

databaseName = 'fulldb_paspeech_speech_tuning_allcells.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
otherDatabaseName = 'fulldb_speech_tuning.h5'
otherDatabaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, otherDatabaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
otherdb = celldatabase.load_hdf(otherDatabaseFullPath)
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']

alpha = 0.05/12
speechResponsive = np.array(celldb.speechMinPvalOnset < alpha)

## -- exclude cells with low firing rates
exclusionCriterion = 3  # sp/s
maxFiringRateSpeechEvoked = np.max([otherdb.maxFiringRate_FT_VOTmin, otherdb.maxFiringRate_FT_VOTmax, otherdb.maxFiringRate_VOT_FTmin, otherdb.maxFiringRate_VOT_FTmax], 0)
exclude_evoked = maxFiringRateSpeechEvoked < exclusionCriterion
exclude_baseline = otherdb.speechFiringRateBaseline < exclusionCriterion
excludeCells = exclude_evoked & exclude_baseline


## -- get min pval b/w conditions
ftCombos = np.array([celldb.ftSelectivityVotMaxPvalOnset, celldb.ftSelectivityVotMinPvalOnset])
votCombos = np.array([celldb.votSelectivityFtMaxPvalOnset, celldb.votSelectivityFtMinPvalOnset])
bestSelectivityFt = ftCombos.min(0)
bestSelectivityVot = votCombos.min(0)

amOnsetSelective = np.array(celldb.amSelectivityPvalOnset < 0.05)
amSustainSelective = np.array(celldb.amSelectivityPvalSustain < 0.05)
amSelective = amOnsetSelective| amSustainSelective
toneSelective = np.array(celldb.toneSelectivityPval < 0.05)
excludeAM = celldb.amFiringRateMaxOnset < exclusionCriterion
excludeTone = celldb.toneFiringRateBest < exclusionCriterion
amSelective[excludeAM] = np.nan
toneSelective[excludeTone] = np.nan
deepLayers = celldb.recordingSiteName.str.contains('layer 5') | celldb.recordingSiteName.str.contains('layer 6')
superLayers = celldb.recordingSiteName.str.contains('layer 1') | celldb.recordingSiteName.str.contains('layer 2') | celldb.recordingSiteName.str.contains('layer 3')
layer4 = celldb.recordingSiteName.str.contains('layer 4')

## -- set selectivityCriterion for mixed selectivity
selectivityCriterion = 0.025 #0.05/2 (because taking the best selectivity across two conditions)
VOTselective = bestSelectivityVot < selectivityCriterion
FTselective = bestSelectivityFt < selectivityCriterion
mixedSelective = VOTselective & FTselective
singleSelective = np.logical_xor(VOTselective, FTselective)
notSelective = ~VOTselective & ~FTselective
y_coords = celldb.y_coord.copy()
z_coords = celldb.z_coord.copy()

FtminPvalbyArea = []
VotminPvalbyArea = []
speechResponsiveByArea = []
amSelectivebyArea = []
toneSelectivebyArea = []
excludeCellsbyArea = []
mixedSelectivebyArea = []
singleSelectivebyArea = []
notSelectivebyArea = []
VOTselectivebyArea = []
FTselectivebyArea = []
deepLayersbyArea = []
superLayersbyArea = []
layer4byArea = []



if AREASBYDEPTH == 0: ## used for original analysis, where defined areas by allen atlas area
    recordingAreaName = celldb.copy()['recordingAreaName']
    recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
    recordingAreaName = np.array(recordingAreaName)
    for indArea, thisArea in enumerate(audCtxAreas):
        FtminPvalbyArea.append(bestSelectivityFt[recordingAreaName == thisArea])
        VotminPvalbyArea.append(bestSelectivityVot[recordingAreaName == thisArea])
        speechResponsiveByArea.append(speechResponsive[recordingAreaName == thisArea])
        amSelectivebyArea.append(amSelective[recordingAreaName == thisArea])
        toneSelectivebyArea.append(toneSelective[recordingAreaName == thisArea])
        excludeCellsbyArea.append(excludeCells[recordingAreaName == thisArea])
        mixedSelectivebyArea.append(mixedSelective[recordingAreaName == thisArea])
        singleSelectivebyArea.append(singleSelective[recordingAreaName == thisArea])
        notSelectivebyArea.append(notSelective[recordingAreaName == thisArea])
        VOTselectivebyArea.append(VOTselective[recordingAreaName == thisArea])
        FTselectivebyArea.append(FTselective[recordingAreaName == thisArea])
        deepLayersbyArea.append(deepLayers[recordingAreaName == thisArea])
        superLayersbyArea.append(superLayers[recordingAreaName == thisArea])
        layer4byArea.append(layer4[recordingAreaName == thisArea])


if AREASBYDEPTH: ## used for test where we defined ACtx area by depth
    print('using depth to define auditory areas')
    isAudArea = celldb.recordingAreaName.str.contains('auditory area')
    # Coords as seen as splits on hist: audD < 105 < audP < 133 < audV
    # Coords as even thirds: audD < 110 < audP < 124 < audV
    '''
    audP_coords = (celldb.y_coord > 106) & (celldb.y_coord <= 124) #& isAudArea #& (celldb.z_coord > 185)
    audD_coords = (celldb.y_coord <= 106) #& isAudArea #& (celldb.z_coord > 185)
    audV_coords = (celldb.y_coord > 124) #& isAudArea #& (celldb.z_coord > 185)
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
        FtminPvalbyArea.append(bestSelectivityFt[audAreas_byCoords[indArea]])
        VotminPvalbyArea.append(bestSelectivityVot[audAreas_byCoords[indArea]])
        speechResponsiveByArea.append(speechResponsive[audAreas_byCoords[indArea]])
        amSelectivebyArea.append(amSelective[audAreas_byCoords[indArea]])
        toneSelectivebyArea.append(toneSelective[audAreas_byCoords[indArea]])
        excludeCellsbyArea.append(excludeCells[audAreas_byCoords[indArea]])
        mixedSelectivebyArea.append(mixedSelective[audAreas_byCoords[indArea]])
        singleSelectivebyArea.append(singleSelective[audAreas_byCoords[indArea]])
        notSelectivebyArea.append(notSelective[audAreas_byCoords[indArea]])
        VOTselectivebyArea.append(VOTselective[audAreas_byCoords[indArea]])
        FTselectivebyArea.append(FTselective[audAreas_byCoords[indArea]])



## -- exclude low spike count cells
for indArea, thisArea in enumerate(audCtxAreas):
    FtminPvalbyArea[indArea] = FtminPvalbyArea[indArea][~excludeCellsbyArea[indArea]]
    VotminPvalbyArea[indArea] = VotminPvalbyArea[indArea][~excludeCellsbyArea[indArea]]
    speechResponsiveByArea[indArea] = speechResponsiveByArea[indArea][~excludeCellsbyArea[indArea]]
    amSelectivebyArea[indArea] = amSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    toneSelectivebyArea[indArea] = toneSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    mixedSelectivebyArea[indArea] = mixedSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    singleSelectivebyArea[indArea] = singleSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    notSelectivebyArea[indArea] = notSelectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    VOTselectivebyArea[indArea] = VOTselectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    FTselectivebyArea[indArea] = FTselectivebyArea[indArea][~excludeCellsbyArea[indArea]]
    if AREASBYDEPTH == 0:
        deepLayersbyArea[indArea] = deepLayersbyArea[indArea][~excludeCellsbyArea[indArea]]
        superLayersbyArea[indArea] = superLayersbyArea[indArea][~excludeCellsbyArea[indArea]]
        layer4byArea[indArea] = layer4byArea[indArea][~excludeCellsbyArea[indArea]]


## -- run stats
# -- test selectivity distribution across cortical areas
oddsratio, pvalFracSelective_AudPvsAudD = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[1])]]))
oddsratio, pvalFracSelective_AudPvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[2])]]))
oddsratio, pvalFracSelective_AudDvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[2]), np.sum(notSelectivebyArea[1])]]))

## -- test mixed selectivity
oddsratio, pvalFracMixed_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[1])]]))
oddsratio, pvalFracMixed_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[2])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[2])]]))
oddsratio, pvalFracMixed_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[2]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[2]), np.sum(singleSelectivebyArea[1])]]))

## -- test each feature selectivity (responsive cells)
oddsratio, pvalFracFtSel_AudPvsAudD = stats.fisher_exact(np.array([[np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])],[np.sum(speechResponsiveByArea[0]) -np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[1]) - np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])]]))
oddsratio, pvalFracFtSel_AudPvsAudV = stats.fisher_exact(np.array([[np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])],[np.sum(speechResponsiveByArea[0]) -np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[2]) - np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])]]))
oddsratio, pvalFracFtSel_AudDvsAudV = stats.fisher_exact(np.array([[np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1]), np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])],[np.sum(speechResponsiveByArea[1]) -np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1]), np.sum(speechResponsiveByArea[2]) - np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])]]))

oddsratio, pvalFracVotSel_AudPvsAudD = stats.fisher_exact(np.array([[np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])],[np.sum(speechResponsiveByArea[0]) -np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[1]) - np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])]]))
oddsratio, pvalFracVotSel_AudPvsAudV = stats.fisher_exact(np.array([[np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])],[np.sum(speechResponsiveByArea[0]) -np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0]), np.sum(speechResponsiveByArea[2]) - np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])]]))
oddsratio, pvalFracVotSel_AudDvsAudV = stats.fisher_exact(np.array([[np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1]), np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])],[np.sum(speechResponsiveByArea[1]) -np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1]), np.sum(speechResponsiveByArea[2]) - np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])]]))

## -- test each feature selectivity (all cells)
oddsratio, pvalFracFtSel_AudPvsAudD_allcells = stats.fisher_exact(np.array([[np.sum(FtminPvalbyArea[0]<0.025), np.sum(FtminPvalbyArea[1]<0.025)],[len(FtminPvalbyArea[0]) -np.sum(FtminPvalbyArea[0]<0.025), len(FtminPvalbyArea[1]) -np.sum(FtminPvalbyArea[1]<0.025)]]))
oddsratio, pvalFracFtSel_AudPvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(FtminPvalbyArea[0]<0.025), np.sum(FtminPvalbyArea[2]<0.025)],[len(FtminPvalbyArea[0]) -np.sum(FtminPvalbyArea[0]<0.025), len(FtminPvalbyArea[2]) -np.sum(FtminPvalbyArea[2]<0.025)]]))
oddsratio, pvalFracFtSel_AudDvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(FtminPvalbyArea[2]<0.025), np.sum(FtminPvalbyArea[1]<0.025)],[len(FtminPvalbyArea[2]) -np.sum(FtminPvalbyArea[2]<0.025), len(FtminPvalbyArea[1]) -np.sum(FtminPvalbyArea[1]<0.025)]]))

oddsratio, pvalFracVotSel_AudPvsAudD_allcells = stats.fisher_exact(np.array([[np.sum(VotminPvalbyArea[0]<0.025), np.sum(VotminPvalbyArea[1]<0.025)],[len(VotminPvalbyArea[0]) -np.sum(VotminPvalbyArea[0]<0.025), len(VotminPvalbyArea[1]) -np.sum(VotminPvalbyArea[1]<0.025)]]))
oddsratio, pvalFracVotSel_AudPvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(VotminPvalbyArea[0]<0.025), np.sum(VotminPvalbyArea[2]<0.025)],[len(VotminPvalbyArea[0]) -np.sum(VotminPvalbyArea[0]<0.025), len(VotminPvalbyArea[2]) -np.sum(VotminPvalbyArea[2]<0.025)]]))
oddsratio, pvalFracVotSel_AudDvsAudV_allcells = stats.fisher_exact(np.array([[np.sum(VotminPvalbyArea[2]<0.025), np.sum(VotminPvalbyArea[1]<0.025)],[len(VotminPvalbyArea[2]) -np.sum(VotminPvalbyArea[2]<0.025), len(VotminPvalbyArea[1]) -np.sum(VotminPvalbyArea[1]<0.025)]]))


if 0: #test by layer
    ## -- test mixed selectivity
    oddsratio, pvalFracMixed_AudPvsAudD_deep = stats.fisher_exact(np.array([[np.sum((mixedSelectivebyArea[0] & deepLayersbyArea[0])), np.sum((mixedSelectivebyArea[1] & deepLayersbyArea[1]))],[np.sum((singleSelectivebyArea[0] & deepLayersbyArea[0])), np.sum((singleSelectivebyArea[1] & deepLayersbyArea[1]))]]))
    oddsratio, pvalFracMixed_AudPvsAudV_deep = stats.fisher_exact(np.array([[np.sum((mixedSelectivebyArea[0] & deepLayersbyArea[0])), np.sum((mixedSelectivebyArea[2] & deepLayersbyArea[2]))],[np.sum((singleSelectivebyArea[0] & deepLayersbyArea[0])), np.sum((singleSelectivebyArea[2] & deepLayersbyArea[2]))]]))
    oddsratio, pvalFracMixed_AudDvsAudV_deep = stats.fisher_exact(np.array([[np.sum((mixedSelectivebyArea[2] & deepLayersbyArea[2])), np.sum((mixedSelectivebyArea[1] & deepLayersbyArea[1]))],[np.sum((singleSelectivebyArea[2] & deepLayersbyArea[2])), np.sum((singleSelectivebyArea[1] & deepLayersbyArea[1]))]]))

    ## -- test each feature selectivity (responsive cells)
    FtRespAudP = np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0] & deepLayersbyArea[0])
    FtRespAudD = np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1] & deepLayersbyArea[1])
    FtRespAudV = np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2] & deepLayersbyArea[2])
    VotRespAudP = np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0] & deepLayersbyArea[0])
    VotRespAudD = np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1] & deepLayersbyArea[1])
    VotRespAudV = np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2] & deepLayersbyArea[2])


    oddsratio, pvalFracFtSel_AudPvsAudD_deep = stats.fisher_exact(np.array([[FtRespAudP, FtRespAudD],[np.sum(speechResponsiveByArea[0] & deepLayersbyArea[0]) - FtRespAudP, np.sum(speechResponsiveByArea[1] & deepLayersbyArea[1]) - FtRespAudD]]))
    oddsratio, pvalFracFtSel_AudPvsAudV_deep = stats.fisher_exact(np.array([[FtRespAudP, FtRespAudV],[np.sum(speechResponsiveByArea[0] & deepLayersbyArea[0]) - FtRespAudP, np.sum(speechResponsiveByArea[2] & deepLayersbyArea[2]) - FtRespAudV]]))
    oddsratio, pvalFracFtSel_AudDvsAudV_deep = stats.fisher_exact(np.array([[FtRespAudD, FtRespAudV],[np.sum(speechResponsiveByArea[1] & deepLayersbyArea[1]) - FtRespAudD, np.sum(speechResponsiveByArea[2] & deepLayersbyArea[2]) - FtRespAudV]]))

    oddsratio, pvalFracVotSel_AudPvsAudD_deep = stats.fisher_exact(np.array([[VotRespAudP, VotRespAudD],[np.sum(speechResponsiveByArea[0] & deepLayersbyArea[0]) - VotRespAudP, np.sum(speechResponsiveByArea[1] & deepLayersbyArea[1]) - VotRespAudD]]))
    oddsratio, pvalFracVotSel_AudPvsAudV_deep = stats.fisher_exact(np.array([[VotRespAudP, VotRespAudV],[np.sum(speechResponsiveByArea[0] & deepLayersbyArea[0]) - VotRespAudP, np.sum(speechResponsiveByArea[2] & deepLayersbyArea[2]) - VotRespAudV]]))
    oddsratio, pvalFracVotSel_AudDvsAudV_deep = stats.fisher_exact(np.array([[VotRespAudV, VotRespAudD],[np.sum(speechResponsiveByArea[2] & deepLayersbyArea[2]) - VotRespAudV, np.sum(speechResponsiveByArea[1] & deepLayersbyArea[1]) - VotRespAudD]]))


    ## -- test each feature selectivity (all cells)
    oddsratio, pvalFracFtSel_AudPvsAudD_allcells_deep = stats.fisher_exact(np.array([[FtRespAudP, FtRespAudD],[(np.sum(deepLayersbyArea[0]) - FtRespAudP), (np.sum(deepLayersbyArea[1]) - FtRespAudD)]]))
    oddsratio, pvalFracFtSel_AudPvsAudV_allcells_deep = stats.fisher_exact(np.array([[FtRespAudP, FtRespAudV],[(np.sum(deepLayersbyArea[0]) - FtRespAudP), (np.sum(deepLayersbyArea[2]) - FtRespAudV)]]))
    oddsratio, pvalFracFtSel_AudDvsAudV_allcells_deep = stats.fisher_exact(np.array([[FtRespAudD, FtRespAudV],[(np.sum(deepLayersbyArea[1]) - FtRespAudD), (np.sum(deepLayersbyArea[2]) - FtRespAudV)]]))

    oddsratio, pvalFracVotSel_AudPvsAudD_allcells_deep = stats.fisher_exact(np.array([[VotRespAudP, VotRespAudD],[(np.sum(deepLayersbyArea[0]) - VotRespAudP), (np.sum(deepLayersbyArea[1]) - VotRespAudD)]]))
    oddsratio, pvalFracVotSel_AudPvsAudV_allcells_deep = stats.fisher_exact(np.array([[VotRespAudP, VotRespAudV],[(np.sum(deepLayersbyArea[0]) - VotRespAudP), (np.sum(deepLayersbyArea[2]) - VotRespAudV)]]))
    oddsratio, pvalFracVotSel_AudDvsAudV_allcells_deep = stats.fisher_exact(np.array([[VotRespAudD, VotRespAudV],[(np.sum(deepLayersbyArea[1]) - VotRespAudD), (np.sum(deepLayersbyArea[2]) - VotRespAudV)]]))

    print('--Stats Summary: LAYER RESULTS--')
    print(f'Frac mixed selective AudP vs AudD:  p = {pvalFracMixed_AudPvsAudD_deep}')
    print(f'Frac mixed selective AudP vs AudV:  p = {pvalFracMixed_AudPvsAudV_deep}')
    print(f'Frac mixed selective AudD vs AudV:  p = {pvalFracMixed_AudDvsAudV_deep}')
    print(f'Frac FT selective AudP vs AudD:  p = {pvalFracFtSel_AudPvsAudD_deep}')
    print(f'Frac FT selective AudP vs AudV:  p = {pvalFracFtSel_AudPvsAudV_deep}')
    print(f'Frac FT selective AudD vs AudV:  p = {pvalFracFtSel_AudDvsAudV_deep}')
    print(f'Frac VOT selective AudP vs AudD:  p = {pvalFracVotSel_AudPvsAudD_deep}')
    print(f'Frac VOT selective AudP vs AudV:  p = {pvalFracVotSel_AudPvsAudV_deep}')
    print(f'Frac VOT selective AudD vs AudV:  p = {pvalFracVotSel_AudDvsAudV_deep}')
    print(f'AudP: n = {np.sum(deepLayersbyArea[0])} layers 5/6, n = {np.sum(superLayersbyArea[0])} layers 1/2/3, n = {np.sum(layer4byArea[0])} layer 4')
    print(f'AudD: n = {np.sum(deepLayersbyArea[1])} layers 5/6, n = {np.sum(superLayersbyArea[1])} layers 1/2/3, n = {np.sum(layer4byArea[1])} layer 4')
    print(f'AudV: n = {np.sum(deepLayersbyArea[2])} layers 5/6, n = {np.sum(superLayersbyArea[2])} layers 1/2/3, n = {np.sum(layer4byArea[2])} layer 4')
    print(f'AudP n speech responsive {np.sum(speechResponsiveByArea[0] & deepLayersbyArea[0])}, n selective: VOT = {VotRespAudP}, FT = {FtRespAudP}')
    print(f'AudD n speech responsive {np.sum(speechResponsiveByArea[1] & deepLayersbyArea[1])}, n selective: VOT = {VotRespAudD}, FT = {FtRespAudD}')
    print(f'AudV n speech responsive {np.sum(speechResponsiveByArea[2] & deepLayersbyArea[2])}, n selective: VOT = {VotRespAudV}, FT = {FtRespAudV}')
    #print(f'Frac Selective AudP = {(np.sum(mixedSelectivebyArea[0]) + np.sum(singleSelectivebyArea[0]))/len(singleSelectivebyArea[0])}')
    #print(f'Frac Selective AudD = {(np.sum(mixedSelectivebyArea[1]) + np.sum(singleSelectivebyArea[1]))/len(singleSelectivebyArea[1])}')
    #print(f'Frac Selective AudV = {(np.sum(mixedSelectivebyArea[2]) + np.sum(singleSelectivebyArea[2]))/len(singleSelectivebyArea[2])}')



'''
## -- test correlation b/w Speech feature selectivity and basic sound selectivity
kstat, pvalAmVot_AudP = stats.kruskal(VotminPvalbyArea[0][amSelectivebyArea[0] == 1],VotminPvalbyArea[0][amSelectivebyArea[0] == 0])
kstat, pvalAmVot_AudD = stats.kruskal(VotminPvalbyArea[1][amSelectivebyArea[1] == 1],VotminPvalbyArea[1][amSelectivebyArea[1] == 0])
kstat, pvalAmVot_AudV = stats.kruskal(VotminPvalbyArea[2][amSelectivebyArea[2] == 1],VotminPvalbyArea[2][amSelectivebyArea[2] == 0])

kstat, pvalToneFt_AudP = stats.kruskal(FtminPvalbyArea[0][toneSelectivebyArea[0]==1], FtminPvalbyArea[0][toneSelectivebyArea[0]==0])
kstat, pvalToneFt_AudD = stats.kruskal(FtminPvalbyArea[1][toneSelectivebyArea[1]==1], FtminPvalbyArea[1][toneSelectivebyArea[1]==0])
kstat, pvalToneFt_AudV = stats.kruskal(FtminPvalbyArea[2][toneSelectivebyArea[2]==1], FtminPvalbyArea[2][toneSelectivebyArea[2]==0])
'''


if STATSUMMARY:
    print('--Stats Summary--')
    '''
    print(f'Group differences SI-FT: p = {np.round(pValKruskalBestFT_allcells,3)}')
    print(f'Group differences SI-VOT: p = {np.round(pValKruskalBestVOT_allcells,3)}')
    if pValKruskalBestVOT < 0.05:
        print(f'SI-VOT AudP vs AudD: p = {np.round(pValmannU_votAudPvsAudD_allcells,3)}')
        print(f'SI-VOT AudP vs AudV: p = {np.round(pValmannU_votAudPvsAudV_allcells,3)}')
        print(f'SI-VOT AudD vs AudV: p = {np.round(pValmannU_votAudDvsAudV_allcells,3)}')
    if pValKruskalBestFT < 0.05:
        print(f'SI-FT AudP vs AudD: p = {np.round(pValmannU_ftAudPvsAudD_allcells,3)}')
        print(f'SI-FT AudP vs AudV: p = {np.round(pValmannU_ftAudPvsAudV_allcells,3)}')
        print(f'SI-FT AudD vs AudV: p = {np.round(pValmannU_ftAudDvsAudV_allcells,3)}')
    '''
    print(f'exclusion criterion = firing rate < {exclusionCriterion} sp/s')
    print(f'n Excluded {audCtxAreas[0]}: {np.sum(excludeCellsbyArea[0])}')
    print(f'n Excluded {audCtxAreas[1]}: {np.sum(excludeCellsbyArea[1])}')
    print(f'n Excluded {audCtxAreas[2]}: {np.sum(excludeCellsbyArea[2])}')
    print(f'Selectivity criterion = {selectivityCriterion}')
    print(f'Frac mixed selective AudP vs AudD:  p = {pvalFracMixed_AudPvsAudD}')
    print(f'Frac mixed selective AudP vs AudV:  p = {pvalFracMixed_AudPvsAudV}')
    print(f'Frac mixed selective AudD vs AudV:  p = {pvalFracMixed_AudDvsAudV}')
    print(f'Frac FT selective AudP vs AudD:  p = {pvalFracFtSel_AudPvsAudD}')
    print(f'Frac FT selective AudP vs AudV:  p = {pvalFracFtSel_AudPvsAudV}')
    print(f'Frac FT selective AudD vs AudV:  p = {pvalFracFtSel_AudDvsAudV}')
    print(f'Frac VOT selective AudP vs AudD:  p = {pvalFracVotSel_AudPvsAudD}')
    print(f'Frac VOT selective AudP vs AudV:  p = {pvalFracVotSel_AudPvsAudV}')
    print(f'Frac VOT selective AudD vs AudV:  p = {pvalFracVotSel_AudDvsAudV}')
    print(f'AudP n speechResponsive {np.sum(speechResponsiveByArea[0])}, n selective: VOT = {np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0])} ({np.round((np.sum((VotminPvalbyArea[0]<0.025) & speechResponsiveByArea[0])) /np.sum(speechResponsiveByArea[0])*100,1)}%), FT = {np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0])} ({np.round((np.sum((FtminPvalbyArea[0]<0.025) & speechResponsiveByArea[0])) /np.sum(speechResponsiveByArea[0])*100,1)}%)')
    print(f'AudD n speechResponsive {np.sum(speechResponsiveByArea[1])}, n selective: VOT = {np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])} ({np.round((np.sum((VotminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])) /np.sum(speechResponsiveByArea[1])*100,1)}%), FT = {np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])} ({np.round((np.sum((FtminPvalbyArea[1]<0.025) & speechResponsiveByArea[1])) /np.sum(speechResponsiveByArea[1])*100,1)}%)')
    print(f'AudV n speechResponsive {np.sum(speechResponsiveByArea[2])}, n selective: VOT = {np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])} ({np.round((np.sum((VotminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])) /np.sum(speechResponsiveByArea[2])*100,1)}%), FT = {np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])} ({np.round((np.sum((FtminPvalbyArea[2]<0.025) & speechResponsiveByArea[2])) /np.sum(speechResponsiveByArea[2])*100,1)}%)')





#np.savez(figDataFullPath, bestSelectivityFt = bestSelectivityFt, bestSelectivityVot = bestSelectivityVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, excludeCells = excludeCells, speechResponsive = speechResponsive, amSelective = amSelective, toneSelective = toneSelective, maxFiringRateSpeechEvoked = maxFiringRateSpeechEvoked, y_coord = celldb.y_coord, z_coord = celldb.z_coord)
#print('saved to ' f'{figDataFullPath}')
