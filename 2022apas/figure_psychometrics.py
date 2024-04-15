"""
Figure comparing psychometrics for each cohort.
"""

import os
import sys
import numpy as np
import scipy
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import studyparams
import studyutils
import figparams
from importlib import reload
reload(figparams)
reload(studyutils)

FIGNAME = 'psychometrics_comparison'
figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)

PANELS = [1, 1, 1, 1] # Plot panel i if PANELS[i]==1

SAVE_FIGURE = 1
outputDir = '/tmp/'
figFilename = 'plots_psychometrics' # Do not include extension
figFormat = 'svg' # 'pdf' or 'svg'
#figSize = [15, 6] # In inches (I'm doubling the size)
figSize = [7.5, 4.5] # In inches (I'm doubling the size)

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.01, 0.31, 0.49, 0.68, 0.84]   # Horiz position for panel labels
labelPosY = [0.95, 0.46]    # Vert position for panel labels

# -- Assigned colors (defined in figparams) --
colorEachCond = [figparams.colors['activeOnly'],
                 figparams.colors['activePassive'],
                 figparams.colors['passiveThenActive']]
eachCondLabel = ['A only', 'A + P', 'P : A']
nCond = len(eachCondLabel)

# -- Plot results --
fig = plt.gcf()
fig.clf()
fig.set_facecolor('w')

gsMain = gridspec.GridSpec(2, 2, width_ratios=[0.25, 0.75], height_ratios=[0.5, 0.5])
gsMain.update(left=0.075, right=0.98, top=0.95, bottom=0.1, wspace=0.3, hspace=0.5)

periods = ['early', 'late']
panelLabels = [['A','B','C','D','E'], ['F','G','H','I','J']]


avgCorrectEachCondEachPeriod = []

for indperiod, period in enumerate(periods):
    # -- Load early psycurves --
    #period = 'early'
    #period = 'late'
    psyCurveData = studyutils.load_stage4(period)
    subjects = list(psyCurveData['subjects'])
    fractionLeft = psyCurveData['fractionLeft']
    possibleFMslopes = psyCurveData['possibleFMslope']
    avgCorrect = psyCurveData['avgCorrect']
    psyCurveParams = psyCurveData['psyCurveParams']
    dayRange = psyCurveData['dayRange']
    dayRangeStr = f'{dayRange[0]+1}-{dayRange[-1]}'  # Start from day 1
    correctEachFMslope = [1,1,1,0,0,0] + [-1,-1,-1,1,1,1]*fractionLeft
    correctExtreme = correctEachFMslope[:,[0,5]].mean(axis=1)
    correctNonExtreme = correctEachFMslope[:,1:5].mean(axis=1)
    correctMidStim = correctEachFMslope[:,2:4].mean(axis=1)

    #avgCorrectEachPeriod.append(avgCorrect)
    #subjectsEachPeriod.append(subjects)
    
    fastSubjectsEachCond, eachCond = studyutils.mice_each_condition('fast')
    assert eachCond==['activeOnly', 'activePassive', 'passiveThenActive']
    miceEachCond = fastSubjectsEachCond

    # --- Test a single cohort --
    #miceEachCond[0] = list(set(miceEachCond[0]) & set(studyparams.MICE_ALL_COH4))

    # -- Exclude mice with not enough S3 sessions (did not pass stage3-4 criteria) --
    '''
    excludeMice = ['pamo019']
    for miceThisCond in miceEachCond:
        for mouse in excludeMice:
            if mouse in miceThisCond:
                miceThisCond.remove(mouse)
    '''

    miceIndEachCond = [[], [], []]
    for indcond, mice in enumerate(miceEachCond):
        for mouse in mice:
            miceIndEachCond[indcond].append(subjects.index(mouse))
    nSubjectsEachCond = [len(miceThisCond) for miceThisCond in miceIndEachCond]
    dataEachCond = [fractionLeft[miceIndEachCond[0],:],
                    fractionLeft[miceIndEachCond[1],:],
                    fractionLeft[miceIndEachCond[2],:]]

    # -- Quantities to compare --
    biasEachMouse = np.empty(len(subjects))
    slopeEachMouse = np.empty(len(subjects))
    for inds, subject in enumerate(subjects):
        biasEachMouse[inds] = extrastats.invpsychfun(0.5, *psyCurveParams[inds])
        slopeEachMouse[inds] = studyutils.get_slope(psyCurveParams[inds])
    avgCorrectEachCond = [avgCorrect[miceIndEachCond[0]],
                          avgCorrect[miceIndEachCond[1]],
                          avgCorrect[miceIndEachCond[2]]]
    avgCorrectEachCondEachPeriod.append(avgCorrectEachCond)
    
    '''
    asympPerfAll = 1 - psyCurveParams[:,2:4].mean(axis=1)
    asympPerfEachCond = [asympPerfAll[miceIndEachCond[0]],
                         asympPerfAll[miceIndEachCond[1]],
                         asympPerfAll[miceIndEachCond[2]]]
    '''
    asympPerfEachCond = [ correctExtreme[miceIndEachCond[0]],
                          correctExtreme[miceIndEachCond[1]],
                          correctExtreme[miceIndEachCond[2]] ]
    
    biasEachCond = [biasEachMouse[miceIndEachCond[0]],
                    biasEachMouse[miceIndEachCond[1]],
                    biasEachMouse[miceIndEachCond[2]]]
    absBiasEachCond = [abs(bb) for bb in biasEachCond] 
    slopeEachCond = [slopeEachMouse[miceIndEachCond[0]],
                     slopeEachMouse[miceIndEachCond[1]],
                     slopeEachMouse[miceIndEachCond[2]]]
    perfNonExtreme = [ correctNonExtreme[miceIndEachCond[0]],
                       correctNonExtreme[miceIndEachCond[1]],
                       correctNonExtreme[miceIndEachCond[2]] ]
    perfMidStim = [ correctMidStim[miceIndEachCond[0]],
                       correctMidStim[miceIndEachCond[1]],
                       correctMidStim[miceIndEachCond[2]] ]
    
    # -- Statistics --
    possibleComp = [[0,1], [0,2], [1,2]] # Possible comparisons
    nComparisons = len(possibleComp)
    pValAvgPerf = np.empty(nComparisons)
    pValAsymp = np.empty(nComparisons)
    pValNonExtreme = np.empty(nComparisons)
    pValMidStim = np.empty(nComparisons)
    pValBias = np.empty(nComparisons)
    pValSlope = np.empty(nComparisons)
    for indcomp, thisComp in enumerate(possibleComp):
        print(f'Comparing: {eachCond[thisComp[0]]} to {eachCond[thisComp[1]]}')
        wstat, pValAvgPerf[indcomp] = stats.ranksums(avgCorrectEachCond[thisComp[0]],
                                                     avgCorrectEachCond[thisComp[1]])
        print(f'p-value (avg perf): {pValAvgPerf[indcomp]:0.04f}')
        wstat, pValAsymp[indcomp] = stats.ranksums(asympPerfEachCond[thisComp[0]],
                                                   asympPerfEachCond[thisComp[1]])
        print(f'p-value (asymp perf): {pValAsymp[indcomp]:0.04f}')
        wstat, pValNonExtreme[indcomp] = stats.ranksums(perfNonExtreme[thisComp[0]],
                                                        perfNonExtreme[thisComp[1]])
        print(f'p-value (non-extreme perf): {pValNonExtreme[indcomp]:0.04f}')
        wstat, pValMidStim[indcomp] = stats.ranksums(perfMidStim[thisComp[0]],
                                                     perfMidStim[thisComp[1]])
        print(f'p-value (mid-stim perf): {pValMidStim[indcomp]:0.04f}')
        wstat, pValBias[indcomp] = stats.ranksums(absBiasEachCond[thisComp[0]],
                                                  absBiasEachCond[thisComp[1]])
        print(f'p-value (abs bias): {pValBias[indcomp]:0.04f}')
        wstat, pValSlope[indcomp] = stats.ranksums(slopeEachCond[thisComp[0]],
                                                   slopeEachCond[thisComp[1]])
        print(f'p-value (slope): {pValSlope[indcomp]:0.04f}')
        print('')

    # -- Panel: Average psychometric early --
    ax0 = plt.subplot(gsMain[indperiod, 0])
    ax0.annotate(panelLabels[indperiod][0], xy=(labelPosX[0],labelPosY[indperiod]),
                 xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
    if PANELS[0]:
        xPad = 0.2 * (possibleFMslopes[-1] - possibleFMslopes[0])
        fitxval = np.linspace(possibleFMslopes[0]-xPad, possibleFMslopes[-1]+xPad, 40)
        lineEachCond = []
        for indcond, cond in enumerate(eachCond):
            fractionLeftEachValue = dataEachCond[indcond].mean(axis=0)
            fractionLeftSEM = dataEachCond[indcond].std(axis=0)/np.sqrt(nSubjectsEachCond[indcond])
            ciLeftEachValue = np.vstack((fractionLeftEachValue+fractionLeftSEM,
                                         fractionLeftEachValue-fractionLeftSEM))

            # -- Fit sigmoidal curve --
            par0 = [0, 0.5, 0, 0]
            bounds = [[-np.inf, 0.08, 0, 0], [np.inf, np.inf, 0.5, 0.5]]
            curveParams, pcov = scipy.optimize.curve_fit(extrastats.psychfun, possibleFMslopes,
                                                         fractionLeftEachValue, p0=par0, bounds=bounds)
            fityval = extrastats.psychfun(fitxval, *curveParams)
            hfit, = plt.plot(fitxval, 100*fityval, '-', linewidth=2, color=colorEachCond[indcond])
            lineEachCond.append(hfit)
            (pline, pcaps, pbars, pdots) = studyutils.plot_psychometric(possibleFMslopes,
                                                                        fractionLeftEachValue,
                                                                        ciLeftEachValue)
            plt.setp(pcaps, color=colorEachCond[indcond])
            plt.setp(pbars, color=colorEachCond[indcond])
            plt.setp(pdots, mfc=colorEachCond[indcond], mec='none', ms=6)
            pline.set_visible(False)
        plt.xlim([-6, 6])
        plt.ylabel('Leftward choice (%)', fontsize=fontSizeLabels)
        plt.xlabel('FM slope (oct/s)', fontsize=fontSizeLabels)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        plt.show()
        #legendLabels = [f'{lab} ({n} mice)' for lab,n in zip(eachCondLabel,nSubjectsEachCond)]
        plt.text(-5, 90, f'Days {dayRangeStr}', color='0.0', ha='left', va='center', fontsize=fontSizeLabels)
        plt.legend(lineEachCond, eachCondLabel, loc='lower right', fontsize=fontSizeLabels, frameon=False)


    # -- Panel: comparison psychometrics early --
    if PANELS[1]:
        markerSizeComp = 4
        medianWidth = 0.3
        gsComp = gsMain[indperiod,1].subgridspec(1, 4, wspace=1.2)
        axCompList = [0,0,0,0]

        axCompList[0] = plt.subplot(gsComp[0, 0])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = 100*avgCorrectEachCond[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 2)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
            plt.ylabel(f'Avg. performance across\nall stimuli (% correct)', fontsize=fontSizeLabels)
        plt.ylim([50, 100])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValAvgPerf[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 96, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValAvgPerf[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 99, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValAvgPerf[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 96, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[0].annotate(panelLabels[indperiod][1], xy=(labelPosX[1],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

        # -- Plot extreme stim (asymptotic) performance --
        axCompList[1] = plt.subplot(gsComp[0, 1])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = 100*asympPerfEachCond[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 1)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp, clip_on=False)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
        plt.ylabel(f'Performance extreme\nstim (% correct)', fontsize=fontSizeLabels)
        plt.ylim([50, 100])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValAsymp[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 96+0, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValAsymp[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 99+0, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValAsymp[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 96+0, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[1].annotate(panelLabels[indperiod][2], xy=(labelPosX[2],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

        # -- Plot performance for non-extreme stim --
        axCompList[2] = plt.subplot(gsComp[0, 2])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = 100*perfNonExtreme[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 1)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp, clip_on=False)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
        plt.ylabel(f'Performance non-extreme\nstim (% correct)', fontsize=fontSizeLabels)
        plt.ylim([50, 100])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValNonExtreme[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 96+0, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValNonExtreme[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 99+0, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValNonExtreme[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 96+0, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[2].annotate(panelLabels[indperiod][3], xy=(labelPosX[3],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

        '''
        # -- Plot performance for mid stim --
        axCompList[3] = plt.subplot(gsComp[0, 3])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = 100*perfMidStim[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 1)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp, clip_on=False)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
        plt.ylabel(f'Performance non-extreme\nstim (% correct)', fontsize=fontSizeLabels)
        plt.ylim([50, 100])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValMidStim[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 96+5, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValMidStim[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 99+5, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValMidStim[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 96+5, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[1].annotate(panelLabels[indperiod][2], xy=(labelPosX[2],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
        '''
        '''
        # -- Plot bias --
        axCompList[2] = plt.subplot(gsComp[0, 2])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = absBiasEachCond[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 0.05)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
        plt.ylabel(f'Abs bias (oct/s)', fontsize=fontSizeLabels)
        #plt.ylim([-0.5, 4.5])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValBias[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 5, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValBias[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 4, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValBias[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 5, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[2].annotate(panelLabels[indperiod][3], xy=(labelPosX[3],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')
        '''
        # -- Plot slope --
        axCompList[3] = plt.subplot(gsComp[0, 3])
        for indcond, cond in enumerate(eachCond):
            dataToPlot = 100*slopeEachCond[indcond]
            medianThisCond = np.median(dataToPlot)
            offset = extraplots.spread_offsets(dataToPlot, 0.1, 2)
            plt.plot(np.tile(indcond, nSubjectsEachCond[indcond])+offset, dataToPlot,
                     'o', mfc='none', mec=colorEachCond[indcond], mew=1, ms=markerSizeComp)
            plt.plot(indcond+medianWidth*np.array([-1,1]), np.tile(medianThisCond,2), lw=3,
                     color=colorEachCond[indcond])
        plt.ylabel(f'Psychometric slope\n(% / (oct/s))', fontsize=fontSizeLabels)
        plt.ylim([0, 40])
        plt.xlim([-0.75, 2.75])
        plt.xticks([0, 1, 2], eachCondLabel, fontsize=fontSizeLabels, rotation=45)
        if pValSlope[0]<0.05:
            hs, hl = extraplots.significance_stars([0, 1], 40, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValSlope[1]<0.05:
            hs, hl = extraplots.significance_stars([0, 2], 38, yLength=0.9, gapFactor=0.2, starSize=6)
            plt.setp(hl, lw=0.75)
        if pValSlope[2]<0.05:
            hs, hl = extraplots.significance_stars([1, 2], 40, yLength=0.9, gapFactor=0.3, starSize=6)
            plt.setp(hl, lw=0.75)
        extraplots.set_ticks_fontsize(plt.gca(), fontSizeTicks)
        extraplots.boxoff(plt.gca())
        axCompList[3].annotate(panelLabels[indperiod][4], xy=(labelPosX[4],labelPosY[indperiod]),
                               xycoords='figure fraction', fontsize=fontSizePanel, fontweight='bold')

        
    plt.show()

avgCorrectEachCondEachPeriod = [np.stack((avgCorrectEachCondEachPeriod[0][indc],
                                          avgCorrectEachCondEachPeriod[1][indc])) for indc in range(nCond)]

print('Comparison of Avg performance between early and late')
for indc in range(nCond):
    (avgCorrectEachCondEachPeriod[indc][0], avgCorrectEachCondEachPeriod[indc][1])
    wtat, pVal = stats.wilcoxon(avgCorrectEachCondEachPeriod[indc][0], avgCorrectEachCondEachPeriod[indc][1])
    print(f'{eachCondLabel[indc]} p={pVal}')

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)

'''
plt.figure(2)
plt.clf()
for indc in range(nCond):
    plt.subplot(1,nCond,indc+1)
    plt.plot(100*avgCorrectEachCondEachPeriod[indc],'-o')
    plt.ylim(50,100)
    plt.ylabel('Average Perf (%)')
    plt.ylabel('Period')
    plt.title(eachCondLabel[indc])
    plt.xticks([0,1], periods)
plt.show()
'''
