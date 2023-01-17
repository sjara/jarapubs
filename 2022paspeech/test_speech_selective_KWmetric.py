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
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)

'''
selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])
'''
alpha = 0.05/12
speechResponsive = np.array(celldb.speechMinPvalOnset < alpha)

## -- exclude cells with low firing rates
exclusionCriterion = 4  # sp/s
maxFiringRateSpeechEvoked = np.max([otherdb.maxFiringRate_FT_VOTmin, otherdb.maxFiringRate_FT_VOTmax, otherdb.maxFiringRate_VOT_FTmin, otherdb.maxFiringRate_VOT_FTmax], 0)
exclude_evoked = maxFiringRateSpeechEvoked < exclusionCriterion
exclude_baseline = otherdb.speechFiringRateBaseline < exclusionCriterion
excludeCells = exclude_evoked & exclude_baseline

## -- select best index b/w conditions
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

## -- set selectivityCriterion for mixed selectivity
selectivityCriterion = 0.025
VOTselective = bestSelectivityVot < selectivityCriterion
FTselective = bestSelectivityFt < selectivityCriterion
mixedSelective = VOTselective & FTselective
singleSelective = np.logical_xor(VOTselective, FTselective)
notSelective = ~VOTselective & ~FTselective

bestFtSIbyArea = []
bestVotSIbyArea = []
speechResponsiveByArea = []
amSelectivebyArea = []
toneSelectivebyArea = []
excludeCellsbyArea = []
mixedSelectivebyArea = []
singleSelectivebyArea = []
notSelectivebyArea = []

for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea.append(bestSelectivityFt[recordingAreaName == thisArea])
    bestVotSIbyArea.append(bestSelectivityVot[recordingAreaName == thisArea])
    speechResponsiveByArea.append(speechResponsive[recordingAreaName == thisArea])
    amSelectivebyArea.append(amSelective[recordingAreaName == thisArea])
    toneSelectivebyArea.append(toneSelective[recordingAreaName == thisArea])
    excludeCellsbyArea.append(excludeCells[recordingAreaName == thisArea])
    mixedSelectivebyArea.append(mixedSelective[recordingAreaName == thisArea])
    singleSelectivebyArea.append(singleSelective[recordingAreaName == thisArea])
    notSelectivebyArea.append(notSelective[recordingAreaName == thisArea])


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


## -- run stats
'''
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
'''
## -- test selectivity distribution across cortical areas
oddsratio, pvalFracSelective_AudPvsAudD = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[1])]]))
oddsratio, pvalFracSelective_AudPvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[0]) + np.sum(mixedSelectivebyArea[0])), (np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2]))],[np.sum(notSelectivebyArea[0]), np.sum(notSelectivebyArea[2])]]))
oddsratio, pvalFracSelective_AudDvsAudV = stats.fisher_exact(np.array([[(np.sum(singleSelectivebyArea[2]) + np.sum(mixedSelectivebyArea[2])), (np.sum(singleSelectivebyArea[1]) + np.sum(mixedSelectivebyArea[1]))],[np.sum(notSelectivebyArea[2]), np.sum(notSelectivebyArea[1])]]))

## -- test mixed selectivity
oddsratio, pvalFracMixed_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[1])]]))
oddsratio, pvalFracMixed_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[0]), np.sum(mixedSelectivebyArea[2])],[np.sum(singleSelectivebyArea[0]), np.sum(singleSelectivebyArea[2])]]))
oddsratio, pvalFracMixed_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(mixedSelectivebyArea[2]), np.sum(mixedSelectivebyArea[1])],[np.sum(singleSelectivebyArea[2]), np.sum(singleSelectivebyArea[1])]]))

## -- test each feature selectivity
oddsratio, pvalFracFtSel_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(bestFtSIbyArea[0]<0.025), np.sum(bestFtSIbyArea[1]<0.025)],[len(bestFtSIbyArea[0]) -np.sum(bestFtSIbyArea[0]<0.025), len(bestFtSIbyArea[1]) -np.sum(bestFtSIbyArea[1]<0.025)]]))
oddsratio, pvalFracFtSel_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(bestFtSIbyArea[0]<0.025), np.sum(bestFtSIbyArea[2]<0.025)],[len(bestFtSIbyArea[0]) -np.sum(bestFtSIbyArea[0]<0.025), len(bestFtSIbyArea[2]) -np.sum(bestFtSIbyArea[2]<0.025)]]))
oddsratio, pvalFracFtSel_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(bestFtSIbyArea[2]<0.025), np.sum(bestFtSIbyArea[1]<0.025)],[len(bestFtSIbyArea[2]) -np.sum(bestFtSIbyArea[2]<0.025), len(bestFtSIbyArea[1]) -np.sum(bestFtSIbyArea[1]<0.025)]]))

oddsratio, pvalFracVotSel_AudPvsAudD = stats.fisher_exact(np.array([[np.sum(bestVotSIbyArea[0]<0.025), np.sum(bestVotSIbyArea[1]<0.025)],[len(bestVotSIbyArea[0]) -np.sum(bestVotSIbyArea[0]<0.025), len(bestVotSIbyArea[1]) -np.sum(bestVotSIbyArea[1]<0.025)]]))
oddsratio, pvalFracVotSel_AudPvsAudV = stats.fisher_exact(np.array([[np.sum(bestVotSIbyArea[0]<0.025), np.sum(bestVotSIbyArea[2]<0.025)],[len(bestVotSIbyArea[0]) -np.sum(bestVotSIbyArea[0]<0.025), len(bestVotSIbyArea[2]) -np.sum(bestVotSIbyArea[2]<0.025)]]))
oddsratio, pvalFracVotSel_AudDvsAudV = stats.fisher_exact(np.array([[np.sum(bestVotSIbyArea[2]<0.025), np.sum(bestVotSIbyArea[1]<0.025)],[len(bestVotSIbyArea[2]) -np.sum(bestVotSIbyArea[2]<0.025), len(bestVotSIbyArea[1]) -np.sum(bestVotSIbyArea[1]<0.025)]]))

'''
## -- test correlation b/w Speech feature selectivity and basic sound selectivity
kstat, pvalAmVot_AudP = stats.kruskal(bestVotSIbyArea[0][amSelectivebyArea[0] == 1],bestVotSIbyArea[0][amSelectivebyArea[0] == 0])
kstat, pvalAmVot_AudD = stats.kruskal(bestVotSIbyArea[1][amSelectivebyArea[1] == 1],bestVotSIbyArea[1][amSelectivebyArea[1] == 0])
kstat, pvalAmVot_AudV = stats.kruskal(bestVotSIbyArea[2][amSelectivebyArea[2] == 1],bestVotSIbyArea[2][amSelectivebyArea[2] == 0])

kstat, pvalToneFt_AudP = stats.kruskal(bestFtSIbyArea[0][toneSelectivebyArea[0]==1], bestFtSIbyArea[0][toneSelectivebyArea[0]==0])
kstat, pvalToneFt_AudD = stats.kruskal(bestFtSIbyArea[1][toneSelectivebyArea[1]==1], bestFtSIbyArea[1][toneSelectivebyArea[1]==0])
kstat, pvalToneFt_AudV = stats.kruskal(bestFtSIbyArea[2][toneSelectivebyArea[2]==1], bestFtSIbyArea[2][toneSelectivebyArea[2]==0])
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
    print(f'Frac selective AudP vs AudD: p = {pvalFracSelective_AudPvsAudD}')
    print(f'Frac selective AudP vs AudV: p = {pvalFracSelective_AudPvsAudV}')
    print(f'Frac selective AudD vs AudV: p = {pvalFracSelective_AudDvsAudV}')
    print(f'Frac mixed selective AudP vs AudD:  p = {pvalFracMixed_AudPvsAudD}')
    print(f'Frac mixed selective AudP vs AudV:  p = {pvalFracMixed_AudPvsAudV}')
    print(f'Frac mixed selective AudD vs AudV:  p = {pvalFracMixed_AudDvsAudV}')
    print(f'Frac FT selective AudP vs AudD:  p = {pvalFracFtSel_AudPvsAudD}')
    print(f'Frac FT selective AudP vs AudV:  p = {pvalFracFtSel_AudPvsAudV}')
    print(f'Frac FT selective AudD vs AudV:  p = {pvalFracFtSel_AudDvsAudV}')
    print(f'Frac VOT selective AudP vs AudD:  p = {pvalFracVotSel_AudPvsAudD}')
    print(f'Frac VOT selective AudP vs AudV:  p = {pvalFracVotSel_AudPvsAudV}')
    print(f'Frac VOT selective AudD vs AudV:  p = {pvalFracVotSel_AudDvsAudV}')
    print(f'Frac Selective AudP = {(np.sum(mixedSelectivebyArea[0]) + np.sum(singleSelectivebyArea[0]))/len(singleSelectivebyArea[0])}')
    print(f'Frac Selective AudD = {(np.sum(mixedSelectivebyArea[1]) + np.sum(singleSelectivebyArea[1]))/len(singleSelectivebyArea[1])}')
    print(f'Frac Selective AudV = {(np.sum(mixedSelectivebyArea[2]) + np.sum(singleSelectivebyArea[2]))/len(singleSelectivebyArea[2])}')




#np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityFt = bestSelectivityFt, bestSelectivityVot = bestSelectivityVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, excludeCells = excludeCells, pValKruskalBestFT = pValKruskalBestFT, pValKruskalBestVOT = pValKruskalBestVOT, speechResponsive = speechResponsive, amSelective = amSelective, toneSelective = toneSelective, maxFiringRateSpeechEvoked = maxFiringRateSpeechEvoked)
#print('saved to ' f'{figDataFullPath}')
