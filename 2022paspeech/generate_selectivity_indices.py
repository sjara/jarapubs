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

databaseName = 'fulldb_speech_tuning.h5'
databaseFullPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME, databaseName)
celldb = celldatabase.load_hdf(databaseFullPath)
#audCtxAreas = ['Primary auditory area', 'Posterior auditory area', 'Dorsal auditory area', 'Ventral auditory area']
audCtxAreas = ['Primary auditory area','Dorsal auditory area', 'Ventral auditory area']
recordingAreaName = celldb.recordingAreaName
recordingAreaName = recordingAreaName.str.replace('Posterior auditory area', 'Dorsal auditory area')
recordingAreaName = np.array(recordingAreaName)


selectivityIndexFT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['maxFiringRate_FT_VOTmin'] + celldb['minFiringRate_FT_VOTmin'])
selectivityIndexFT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['maxFiringRate_FT_VOTmax'] + celldb['minFiringRate_FT_VOTmax'])
selectivityIndexVOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['maxFiringRate_VOT_FTmin'] + celldb['minFiringRate_VOT_FTmin'])
selectivityIndexVOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/(celldb['maxFiringRate_VOT_FTmax'] + celldb['minFiringRate_VOT_FTmax'])

selectivityIndex2FT_VOTmin = (celldb['maxFiringRate_FT_VOTmin'] - celldb['minFiringRate_FT_VOTmin'])/(celldb['minFiringRate_FT_VOTmin'])
selectivityIndex2FT_VOTmax = (celldb['maxFiringRate_FT_VOTmax'] - celldb['minFiringRate_FT_VOTmax'])/(celldb['minFiringRate_FT_VOTmax'])
selectivityIndex2VOT_FTmin = (celldb['maxFiringRate_VOT_FTmin'] - celldb['minFiringRate_VOT_FTmin'])/(celldb['minFiringRate_VOT_FTmin'])
selectivityIndex2VOT_FTmax = (celldb['maxFiringRate_VOT_FTmax'] - celldb['minFiringRate_VOT_FTmax'])/( celldb['minFiringRate_VOT_FTmax'])

alpha = 0.05/12
speechResponsive = np.array(celldb.speechMinPvalOnset < alpha)

# -- exclude cells with low firing rates
exclusionCriterion = 3
exclude_Ft_VotMin = celldb.maxFiringRate_FT_VOTmin < exclusionCriterion
exclude_Ft_VotMax = celldb.maxFiringRate_FT_VOTmax < exclusionCriterion
exclude_Vot_FtMin = celldb.maxFiringRate_FT_VOTmin < exclusionCriterion
exclude_Vot_FtMax = celldb.maxFiringRate_FT_VOTmax < exclusionCriterion

selectivityIndexFT_VOTmin[exclude_Ft_VotMin] = 0
selectivityIndexFT_VOTmax[exclude_Ft_VotMax] = 0
selectivityIndexVOT_FTmin[exclude_Vot_FtMin] = 0
selectivityIndexVOT_FTmax[exclude_Vot_FtMax] = 0

# -- select best index b/w conditions
ftCombos = np.array([selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax])
votCombos = np.array([selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax])
bestSelectivityIndexFt = ftCombos.max(0)
bestSelectivityIndexVot = votCombos.max(0)

# -- run stats
bestFtSIbyArea = []
bestVotSIbyArea = []
speechResponsiveByArea = []

for indArea, thisArea in enumerate(audCtxAreas):
    bestFtSIbyArea.append(bestSelectivityIndexFt[recordingAreaName == thisArea])
    bestVotSIbyArea.append(bestSelectivityIndexVot[recordingAreaName == thisArea])
    speechResponsiveByArea.append(speechResponsive[recordingAreaName == thisArea])

## responsive cells
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

## all cells
kstat, pValKruskalBestVOT_allcells = stats.kruskal(bestVotSIbyArea[0], bestVotSIbyArea[1], bestVotSIbyArea[2], nan_policy = 'omit')
kstat, pValKruskalBestFT_allcells = stats.kruskal(bestFtSIbyArea[0], bestFtSIbyArea[1], bestFtSIbyArea[2], nan_policy = 'omit')

# -- if group difference, test individual compariosns:
if pValKruskalBestVOT_allcells < 0.05:
    ustat, pValmannU_votAudPvsAudD_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[1])
    ustat, pValmannU_votAudPvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[0], bestVotSIbyArea[2])
    ustat, pValmannU_votAudDvsAudV_allcells = stats.mannwhitneyu(bestVotSIbyArea[1], bestVotSIbyArea[2])

if pValKruskalBestFT_allcells < 0.05:
    ustat, pValmannU_ftAudPvsAudD_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[1])
    ustat, pValmannU_ftAudPvsAudV_allcells = stats.mannwhitneyu(bestFtSIbyArea[0], bestFtSIbyArea[2])
    ustat, pValmannU_ftAudDvsAudV_allcells = stats.mannwhitneyu(bestGtSIbyArea[1], bestFtSIbyArea[2])

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
    print(f'n Excluded {audCtxAreas[0]}: {np.sum(bestFtSIbyArea[0]==0)}')
    print(f'n Excluded {audCtxAreas[1]}: {np.sum(bestFtSIbyArea[1]==0)}')
    print(f'n Excluded {audCtxAreas[2]}: {np.sum(bestFtSIbyArea[2]==0)}')


np.savez(figDataFullPath, selectivityIndexFT_VOTmin = selectivityIndexFT_VOTmin, selectivityIndexFT_VOTmax = selectivityIndexFT_VOTmax, selectivityIndexVOT_FTmin = selectivityIndexVOT_FTmin, selectivityIndexVOT_FTmax = selectivityIndexVOT_FTmax, bestSelectivityIndexFt = bestSelectivityIndexFt, bestSelectivityIndexVot = bestSelectivityIndexVot, audCtxAreas = audCtxAreas, recordingAreaName = recordingAreaName, exclusionCriterion = exclusionCriterion, pValKruskalBestFT = pValKruskalBestFT, pValKruskalBestVOT = pValKruskalBestVOT, speechResponsive = speechResponsive)
print('saved to ' f'{figDataFullPath}')
