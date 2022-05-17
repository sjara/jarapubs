"""
Compare D1 vs nD1.
"""

import os
import sys
import studyparams
import numpy as np
import pandas as pd
from jaratoolbox import celldatabase
from jaratoolbox import settings
from jaratoolbox import spikesanalysis
from jaratoolbox import ephyscore
from jaratoolbox import behavioranalysis
from scipy import stats
from scipy import signal

import matplotlib.pyplot as plt
from jaratoolbox import extraplots

import studyutils

from importlib import reload
reload(extraplots)
reload(studyutils)


figuresDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME)
dbPath = os.path.join(figuresDataDir, 'astrpi_am_tuning.h5')
celldb = celldatabase.load_hdf(dbPath)
nCells = len(celldb)

restrictND1 = 1  # True: use only ND1 from tetrodes with D1
toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1 = \
    studyutils.select_cells(celldb, restrictND1=restrictND1)

        
cellsWithTone = ~np.isnan(celldb['toneMinPval'])
cellsWithAM = ~np.isnan(celldb['amMinPvalOnset'])
cellsWithSoundSessions = cellsWithTone | cellsWithAM

indD1tone = indD1 & cellsWithTone
indND1tone = indND1 & cellsWithTone
indD1am = indD1 & cellsWithAM
indND1am = indND1 & cellsWithAM
indD1sound = indD1tone | indD1am
indND1sound = indND1tone | indND1am

nCellsWithSoundSessions = np.sum(cellsWithSoundSessions)
nCellsClassifiedWithSoundSession = np.sum(indD1sound | indND1sound)
nCellsD1sound = np.sum(indD1sound)
nCellsND1sound = np.sum(indND1sound)
nCellsD1tone = indD1tone.sum()
nCellsND1tone = indND1tone.sum()
nCellsD1am = indD1am.sum()
nCellsND1am = indND1am.sum()


laserBase = celldb.laserBaseRate200
laserResp = celldb.laserRespRate50
laserRespIndex = (laserResp-laserBase)/(laserResp+laserBase)
laserRespDelta = (laserResp-laserBase)

nSoundResponsiveD1 = np.sum( (toneResponsive & indD1tone) |
                             (amOnsetResponsive & indD1am) |
                             (amSustainResponsive & indD1am) )
nSoundResponsiveND1 = np.sum( (toneResponsive & indND1tone) |
                             (amOnsetResponsive & indND1am) |
                             (amSustainResponsive & indND1am) )
nToneResponsiveD1 = (toneResponsive & indD1tone).sum()
nToneResponsiveND1 = (toneResponsive & indND1tone).sum()
nAMonsetResponsiveD1 = (amOnsetResponsive & indD1am).sum()
nAMonsetResponsiveND1 = (amOnsetResponsive & indND1am).sum()
nAMsustainResponsiveD1 = (amSustainResponsive & indD1am).sum()
nAMsustainResponsiveND1 = (amSustainResponsive & indND1am).sum()

soundResponsiveTable = [ [nSoundResponsiveD1, nSoundResponsiveND1],
                        [nCellsD1tone-nSoundResponsiveD1, nCellsND1tone-nSoundResponsiveND1] ]
oddsratio, soundPvalFisher = stats.fisher_exact(soundResponsiveTable, alternative='two-sided')

toneResponsiveTable = [ [nToneResponsiveD1, nToneResponsiveND1],
                        [nCellsD1tone-nToneResponsiveD1, nCellsND1tone-nToneResponsiveND1] ]
oddsratio, tonePvalFisher = stats.fisher_exact(toneResponsiveTable, alternative='two-sided')

amOnsetResponsiveTable = [ [nAMonsetResponsiveD1, nAMonsetResponsiveND1],
                           [nCellsD1am-nAMonsetResponsiveD1, nCellsND1am-nAMonsetResponsiveND1] ]
oddsratio, amOnsetPvalFisher = stats.fisher_exact(amOnsetResponsiveTable, alternative='two-sided')
amSustainResponsiveTable = [ [nAMsustainResponsiveD1, nAMsustainResponsiveND1],
                           [nCellsD1am-nAMsustainResponsiveD1, nCellsND1am-nAMsustainResponsiveND1] ]
oddsratio, amSustainPvalFisher = stats.fisher_exact(amSustainResponsiveTable, alternative='two-sided')

print(f'Total cells with sound sessions: {nCellsWithSoundSessions}')
print(f'{cellsWithTone.sum()} (tones) and {cellsWithAM.sum()} (AM)')
print()

print(f'Total D1 or ND1 with sound sessions: {nCellsClassifiedWithSoundSession}')
print(f'D1: {nCellsD1sound} and ND1: {nCellsND1sound}')
print()

print(f'Sound responsive: \t' +
      f'D1: {nSoundResponsiveD1}/{nCellsD1sound} ({nSoundResponsiveD1/nCellsD1sound:0.2%})  vs  ' +
      f'nD1: {nSoundResponsiveND1}/{nCellsND1sound} ({nSoundResponsiveND1/nCellsND1sound:0.2%})' +
      f'\t p = {soundPvalFisher:0.4f}')
print()

print(f'N cells (tone): D1:{nCellsD1tone}  vs  nD1:{nCellsND1tone}')

print(f'Tone responsive: \t' +
      f'D1: {nToneResponsiveD1} ({nToneResponsiveD1/nCellsD1tone:0.2%})  vs  ' +
      f'nD1: {nToneResponsiveND1} ({nToneResponsiveND1/nCellsND1tone:0.2%})' +
      f'\t p = {tonePvalFisher:0.4f}')

print()
print(f'N cells (AM) : D1:{nCellsD1am}  vs  nD1:{nCellsND1am}')
print(f'AM Onset responsive: \t' +
      f'D1: {nAMonsetResponsiveD1} ({nAMonsetResponsiveD1/nCellsD1am:0.2%})  vs  ' +
      f'nD1: {nAMonsetResponsiveND1} ({nAMonsetResponsiveND1/nCellsND1am:0.2%})' +
      f'\t p = {amOnsetPvalFisher:0.4f}')

print(f'AM Sustain responsive: \t' +
      f'D1: {nAMsustainResponsiveD1} ({nAMsustainResponsiveD1/nCellsD1am:0.2%})  vs  ' +
      f'nD1: {nAMsustainResponsiveND1} ({nAMsustainResponsiveND1/nCellsND1am:0.2%})' +
      f'\t p = {amSustainPvalFisher:0.4f}')
print()


# -- Tone response index --
toneRespIndex = ( (celldb.toneFiringRateBest-celldb.toneFiringRateBaseline) /
                  (celldb.toneFiringRateBest+celldb.toneFiringRateBaseline) )

if 1:
    print('Using all cells (tone responsive and not-responsive).')
    toneRespIndexD1 = toneRespIndex[indD1]
    toneRespIndexND1 = toneRespIndex[indND1]
else:
    print('Using only TONE RESPONSIVE cells')
    toneRespIndexD1 = toneRespIndex[indD1 & toneResponsive]
    toneRespIndexND1 = toneRespIndex[indND1 & toneResponsive]

toneRespIndexD1pos = toneRespIndexD1[toneRespIndexD1>0]
toneRespIndexND1pos = toneRespIndexND1[toneRespIndexND1>0]
toneRespIndexD1neg = toneRespIndexD1[toneRespIndexD1<0]
toneRespIndexND1neg = toneRespIndexND1[toneRespIndexND1<0]

uval, pValPos = stats.mannwhitneyu(toneRespIndexD1pos, toneRespIndexND1pos, alternative='two-sided')
print(f'Tone index (pos): D1: {toneRespIndexD1pos.mean():0.4f} vs ' +
      f'ND1: {toneRespIndexND1pos.mean():0.4f}   ' +
      f'p = {pValPos:0.6f}')
uval, pValNeg = stats.mannwhitneyu(toneRespIndexD1neg, toneRespIndexND1neg, alternative='two-sided')
print(f'Tone index (neg): D1: {toneRespIndexD1neg.mean():0.4f} vs ' +
      f'ND1: {toneRespIndexND1neg.mean():0.4f}   ' +
      f'p = {pValNeg:0.6f}')
print()


# -- AM Onset response index --
amOnsetRespIndex = ( (celldb.amFiringRateBestOnset-celldb.amFiringRateBaseline) /
                     (celldb.amFiringRateBestOnset+celldb.amFiringRateBaseline) )

if 1:
    print('Using all cells (AM onset responsive and not-responsive).')
    amOnsetRespIndexD1 = amOnsetRespIndex[indD1]
    amOnsetRespIndexND1 = amOnsetRespIndex[indND1]
else:
    print('Using only AM ONSET RESPONSIVE cells')
    amOnsetRespIndexD1 = amOnsetRespIndex[indD1 & amOnsetResponsive]
    amOnsetRespIndexND1 = amOnsetRespIndex[indND1 & amOnsetResponsive]

amOnsetRespIndexD1pos = amOnsetRespIndexD1[amOnsetRespIndexD1>0]
amOnsetRespIndexND1pos = amOnsetRespIndexND1[amOnsetRespIndexND1>0]
amOnsetRespIndexD1neg = amOnsetRespIndexD1[amOnsetRespIndexD1<0]
amOnsetRespIndexND1neg = amOnsetRespIndexND1[amOnsetRespIndexND1<0]

uval, pValPos = stats.mannwhitneyu(amOnsetRespIndexD1pos, amOnsetRespIndexND1pos, alternative='two-sided')
print(f'AmOnset index (pos): D1: {amOnsetRespIndexD1pos.mean():0.4f} vs ' +
      f'ND1: {amOnsetRespIndexND1pos.mean():0.4f}   ' +
      f'p = {pValPos:0.6f}')
uval, pValNeg = stats.mannwhitneyu(amOnsetRespIndexD1neg, amOnsetRespIndexND1neg, alternative='two-sided')
print(f'AmOnset index (neg): D1: {amOnsetRespIndexD1neg.mean():0.4f} vs ' +
      f'ND1: {amOnsetRespIndexND1neg.mean():0.4f}   ' +
      f'p = {pValNeg:0.6f}')
print()


# -- AM Sustain response index --
amSustainRespIndex = ( (celldb.amFiringRateBestSustain-celldb.amFiringRateBaseline) /
                     (celldb.amFiringRateBestSustain+celldb.amFiringRateBaseline) )

if 1:
    print('Using all cells (AM sustain responsive and not-responsive).')
    amSustainRespIndexD1 = amSustainRespIndex[indD1]
    amSustainRespIndexND1 = amSustainRespIndex[indND1]
else:
    print('Using only AM SUSTAIN RESPONSIVE cells')
    amSustainRespIndexD1 = amSustainRespIndex[indD1 & amSustainResponsive]
    amSustainRespIndexND1 = amSustainRespIndex[indND1 & amSustainResponsive]

amSustainRespIndexD1pos = amSustainRespIndexD1[amSustainRespIndexD1>0]
amSustainRespIndexND1pos = amSustainRespIndexND1[amSustainRespIndexND1>0]
amSustainRespIndexD1neg = amSustainRespIndexD1[amSustainRespIndexD1<0]
amSustainRespIndexND1neg = amSustainRespIndexND1[amSustainRespIndexND1<0]

uval, pValPos = stats.mannwhitneyu(amSustainRespIndexD1pos, amSustainRespIndexND1pos, alternative='two-sided')
print(f'AmSustain index (pos): D1: {amSustainRespIndexD1pos.mean():0.4f} vs ' +
      f'ND1: {amSustainRespIndexND1pos.mean():0.4f}   ' +
      f'p = {pValPos:0.6f}')
uval, pValNeg = stats.mannwhitneyu(amSustainRespIndexD1neg, amSustainRespIndexND1neg, alternative='two-sided')
print(f'AmSustain index (neg): D1: {amSustainRespIndexD1neg.mean():0.4f} vs ' +
      f'ND1: {amSustainRespIndexND1neg.mean():0.4f}   ' +
      f'p = {pValNeg:0.6f}')
print()


def line_histogram(data, edges, **kwargs):
    hh, ee = np.histogram(data, edges)
    hline = plt.step(np.r_[ee, ee[-1]], np.r_[0,hh,0],'-', where='pre', **kwargs)
    return hline


plt.clf()
# -- Plot AM sustain response index --
if 1: 
    plt.subplot(3,1,3)
    binEdges = np.linspace(-1, 1, 32)
    line_histogram(amSustainRespIndexD1, binEdges, lw=2, color='b')
    line_histogram(amSustainRespIndexND1, binEdges, lw=2, color='r')
    #plt.hist(amSustainRespIndexD1, binEdges)
    #plt.hist(amSustainRespIndexND1, binEdges)
    plt.title('AM sustained')
    plt.show()
# -- Plot AM onset response index --
if 1: 
    plt.subplot(3,1,2)
    binEdges = np.linspace(-1, 1, 32)
    line_histogram(amOnsetRespIndexD1, binEdges, lw=2, color='b')
    line_histogram(amOnsetRespIndexND1, binEdges, lw=2, color='r')
    #plt.hist(amOnsetRespIndexD1, binEdges)
    #plt.hist(amOnsetRespIndexND1, binEdges)
    plt.title('AM onset')
    plt.show()
# -- Plot tone response index --
if 1:
    plt.subplot(3,1,1)
    binEdges = np.linspace(-1, 1, 32)
    line_histogram(toneRespIndexD1, binEdges, lw=2, color='b')
    line_histogram(toneRespIndexND1, binEdges, lw=2, color='r')
    #plt.hist(toneRespIndexD1, binEdges)
    #plt.hist(toneRespIndexND1, binEdges)
    plt.title('Tones')
    plt.show()
# -- Plot laser response index --
if 0: 
    binEdges = np.linspace(-1, 1, 40)
    plt.hist(laserRespIndex, binEdges)
    plt.hist(laserRespIndex[indD1], binEdges)
    plt.show()
# -- Plot laser response delta --
if 0: 
    binEdges = np.linspace(-50, 100, 200)
    plt.hist(laserRespDelta, binEdges)
    plt.hist(laserRespDelta[indD1], binEdges)
    plt.show()
# -- Plot spike-shapes --
if 0:
    spikeShapes = np.array(list(celldb.spikeShape))

    plt.subplot(2,1,1)
    plt.plot(spikeShapes[indD1,:].T)
    plt.subplot(2,1,2)
    plt.plot(spikeShapes[indND1,:].T)
# -- Plot average spike-shapes --
if 0:
    plt.plot(spikeShapes[indD1,:].mean(axis=0))
    plt.plot(spikeShapes[indND1,:].mean(axis=0))
plt.show()
'''
'''

# -- Show all recording sites --
if 0:
    from jaratoolbox import histologyanalysis as ha
    aa = ha.AllenAverageCoronalAtlas()
    aa.add_points_from_db(celldb)
    aa.show_all_sites()
