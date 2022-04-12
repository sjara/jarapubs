"""
Based on ~/src/jaratest/santiago/talks/201910_tenure_talk/figure_muscimol_psychometric_summary_amod.py
"""

import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats

#dataDir = '/home/nick/data/dissertation_amod/'
dataDir = './'

outputDir = '/tmp/'
figFilename = 'amod_task' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
#figSize = [5, 2] # In inches  2021 submission
figSize = [4, 2] # In inches   2022 submission


# amod002 and 003 ran with chords and AM sounds
subject = 'amod002'
#subject = 'amod003'
soundTypes = ['chords', 'amp_mod']

SHOW_MUSCIMOL = 0


# amod004 ran with tones and AM sounds
# subject = 'amod004'
# soundTypes = ['tones', 'amp_mod']

plt.clf()
gs = gridspec.GridSpec(1,2,width_ratios=[1, 1], height_ratios=[1])
gs.update(top=0.85, left=0.2, right=0.95, bottom=0.25, wspace=0.5, hspace=0.3)

for indSt, stype in enumerate(soundTypes):

    musFilename = 'task_psychometric_{}_{}_muscimol.npz'.format(subject, stype)
    musFullPath = os.path.join(dataDir,musFilename)
    musData = np.load(musFullPath)

    salFilename = 'task_psychometric_{}_{}_saline.npz'.format(subject, stype)
    salFullPath = os.path.join(dataDir,salFilename)
    salData = np.load(salFullPath)

    if SHOW_MUSCIMOL:
        dataToPlot = [musData, salData]
        curveColors = ['r', 'k']
    else:
        dataToPlot = [salData]
        curveColors = ['k']
    plotHandles = []

    fontSizeLabels = 12
    fontSizeTicks = 12
    fontSizePanel = 16

    #ax1 = plt.subplot(1, 2, indSt+1)
    ax1 = plt.subplot(gs[0,indSt])

    for indCond, condData in enumerate(dataToPlot):

        color = curveColors[indCond]

        logPossibleValues = condData['logPossibleValues']
        estimate = condData['estimate']
        ciHitsEachValue = condData['ciHitsEachValue']
        fractionHitsEachValue = condData['fractionHitsEachValue']
        possibleValues = condData['possibleValues']

        xRange = logPossibleValues[-1]-logPossibleValues[1]

        if indSt==0:
            fitxvals = np.linspace(logPossibleValues[0]-0.1*xRange,logPossibleValues[-1]+0.1*xRange,40)
        else:
            edgefactor = 0.2
            fitxvals = np.linspace(logPossibleValues[0]-edgefactor*xRange,logPossibleValues[-1]+edgefactor*xRange,40)
        fityvals = extrastats.psychfun(fitxvals, *estimate)

        upperWhisker = ciHitsEachValue[1,:]-fractionHitsEachValue
        lowerWhisker = fractionHitsEachValue-ciHitsEachValue[0,:]

        # ax1.hold(True)
        (pline, pcaps, pbars) = ax1.errorbar(logPossibleValues,
                                                100*fractionHitsEachValue,
                                                yerr = [100*lowerWhisker, 100*upperWhisker],
                                                ecolor=color, clip_on=False, lw=0) #fmt=None, 

        pdots = ax1.plot(logPossibleValues, 100*fractionHitsEachValue, 'o', ms=6, mec='None', mfc=color,
                         clip_on=False)

        #ax1.set_xticks(logPossibleValues)
        #freqLabels = ['{:.03}'.format(x) for x in possibleValues/1000.0]
        #ax1.set_xticklabels(freqLabels)
        #ax1.set_xlabel('Frequency (kHz)', fontsize=fontSizeLabels)

        pfit, = ax1.plot(fitxvals, 100*fityvals, color=color, lw=2, clip_on=False)
        plotHandles.append(pfit)

    # ax1.annotate('B', xy=(labelPosX[0],labelPosY[0]), xycoords='axes fraction',
    #                 fontsize=fontSizePanel, fontweight='bold')

    extraplots.boxoff(ax1)

    xticks = ax1.get_xticks()
    newXtickLabels = np.logspace(logPossibleValues[0], logPossibleValues[-1], 3, base=2)
    # 1/0
    ax1.set_xticks(np.log2(np.array(newXtickLabels)))
    if not stype=='amp_mod':
        #ax1.set_xticklabels(['{:.3}'.format(x/1000.0) for x in newXtickLabels])
        ax1.set_xticklabels(['{:d}'.format(int(np.round(x/1000.0))) for x in newXtickLabels])
        ax1.set_xlabel('Frequency (kHz)', fontsize=fontSizeLabels)
        ax1.set_title(stype.capitalize())
    else:
        #ax1.set_xticklabels(['{:.3}'.format(x) for x in newXtickLabels])
        ax1.set_xticklabels(['{:d}'.format(int(x)) for x in newXtickLabels])
        ax1.set_xlabel('AM Rate (Hz)', fontsize=fontSizeLabels)
        ax1.set_title("AM")
        #ax1.set_xlim(np.log2([6, 48]))
        ax1.set_xlim(np.log2([5.3, 48]))

    # ax1.set_xlim([fitxvals[0],fitxvals[-1]])

    ax1.set_ylim([0, 100])
    ax1.set_ylabel('Rightward trials (%)', fontsize=fontSizeLabels)
    extraplots.set_ticks_fontsize(plt.gca(),fontSizeTicks)
    ax1.set_yticks([0, 50, 100])

    if SHOW_MUSCIMOL:
        leg = ax1.legend([plotHandles[1],plotHandles[0]], ['Saline','Muscimol'], loc='upper left', frameon=False,
                         labelspacing=0.1, handlelength=1.5, handletextpad=0.2, borderaxespad=0.1, fontsize=12)

plt.show()

#plt.savefig('/home/nick/data/dissertation_amod/figure_{}_muscimol_psychometric_summary.png'.format(subject))

extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
