"""
Generic methods and classes used throughout.
"""

import os
import datetime
import numpy as np
import pandas as pd
import studyparams
from scipy import stats
import matplotlib.pyplot as plt
from jaratoolbox import settings
from jaratoolbox import behavioranalysis
from jaratoolbox import extrastats


def load_sessions_dataframe(subject):
    """
    Load CSV file containing sessions info.
    """
    sessionsDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, 'sessions')
    sessionsFile = os.path.join(sessionsDir, f'{subject}.csv')
    dframe = pd.read_csv(sessionsFile, index_col=0)
    return dframe
    
def get_sessions(subject, stage):
    """
    Return list of sessions for a given stage.
    """
    dframe = load_sessions_dataframe(subject)
    return list(dframe[dframe.stage==stage].session)

def get_antibias(subject, stage):
    """
    Return list of sessions for a given stage.
    """
    dframe = load_sessions_dataframe(subject)
    return list(dframe[(dframe.stage==stage)&dframe.antibias].session)

def rigs_each_subject():
    """
    Show the rigs used for each subject.
    """
    subjects = studyparams.MICE_ALL
    for subject in subjects:
        dframe = load_sessions_dataframe(subject)
        rigs = dframe.rig.unique()
        print(f'{subject}: {rigs}')

def load_stage3(ndays=26, excludeAntibias=False, cohorts=[2,3,4]):
    """
    Load data for stage 3 from specified cohorts into a single dataframe.
    """
    FIGNAME = 'learning_curve_stage3'
    figDataFile = 'fraction_correct_stage3_coh{0}.csv'
    figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
    figDataFullPath = os.path.join(figDataDir, figDataFile)
    #cohorts = [2, 3, 4]
    dframe = pd.DataFrame()
    for indcoh, cohort in enumerate(cohorts):
        dframeOneCohort = pd.read_csv(figDataFullPath.format(cohort), index_col=0)
        # -- Exclude sessions with likely hardware issues --
        for subject, sessions in studyparams.UNRELIABLE_SESSIONS.items():
            for sessionISO in sessions:
                session = sessionISO.replace('-','')+'a'
                if session in dframeOneCohort:
                    dframeOneCohort[session][subject] = np.nan
        if excludeAntibias:
            for subject in dframeOneCohort.index:
                antibiasSessions = get_antibias(subject, 3)
                for abSession in antibiasSessions:
                    dframeOneCohort[abSession][subject] = np.nan
        dframeOneCohort = dframeOneCohort.iloc[:, :ndays]
        dframeOneCohort.columns = np.arange(len(dframeOneCohort.columns))
        dframe = pd.concat((dframe, dframeOneCohort))
    return dframe

def load_stage4(period='early'):
    """
    Load data for stage 4 for all cohorts.
    """
    FIGNAME = 'psychometrics_comparison'
    figDataFile = 'psychometric_{0}_coh{1}.npz'
    figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
    figDataFullPath = os.path.join(figDataDir, figDataFile)
    cohorts = [2, 3, 4]
    # The following parameters come from generate_psychometric_each_mouse.py
    psyCurveDict = {'subjects':[], 'dayRange':None, 'avgCorrect':np.array([]),
                    'possibleFMslope':None, 'fractionLeft':np.array([]),
                    'ciFractionLeft':np.array([]), 'psyCurveParams':np.array([])}
    for indcoh, cohort in enumerate(cohorts):
        datafile = figDataFullPath.format(period,cohort)
        if os.path.isfile(datafile):
            dataOneCohort = np.load(datafile)
        else:
            continue
        if indcoh==0:
            psyCurveDict['subjects'] = dataOneCohort['subjects']
            psyCurveDict['dayRange'] = dataOneCohort['dayRange']
            psyCurveDict['possibleFMslope'] = dataOneCohort['possibleFMslope']
            psyCurveDict['avgCorrect'] = dataOneCohort['avgCorrect']
            psyCurveDict['fractionLeft'] = dataOneCohort['fractionLeft']
            psyCurveDict['ciFractionLeft'] = dataOneCohort['ciFractionLeft']
            psyCurveDict['psyCurveParams'] = dataOneCohort['psyCurveParams']
        else:
            psyCurveDict['subjects'] = np.append(psyCurveDict['subjects'],
                                                   dataOneCohort['subjects'])
            psyCurveDict['avgCorrect'] = np.append(psyCurveDict['avgCorrect'],
                                                   dataOneCohort['avgCorrect'])
            psyCurveDict['fractionLeft'] = np.append(psyCurveDict['fractionLeft'],
                                                     dataOneCohort['fractionLeft'], axis=0)
            psyCurveDict['ciFractionLeft'] = np.append(psyCurveDict['ciFractionLeft'],
                                                       dataOneCohort['ciFractionLeft'], axis=0)
            psyCurveDict['psyCurveParams'] = np.append(psyCurveDict['psyCurveParams'],
                                                       dataOneCohort['psyCurveParams'], axis=0)
    return psyCurveDict
        
def fit_learning_curves(dframe):
    """
    Find linear fit of learning curves.

    Args:
        dframe (pd.DataFrame): pandas dataframe containing fraction correct for each session.
    Returns:
        dframeFit (pd.DataFrame): pandas dataframe with columns related to linear fit
            'slope', 'intercept', 'perfAt21d', 'daysTo70percent'.                       
    """
    nSubjects = len(dframe)
    slopes = np.empty(nSubjects)
    intercepts = np.empty(nSubjects)

    for inds in range(nSubjects):
        perfThisSubject = dframe.iloc[inds].to_numpy()
        # -- Assumes contiguous dates and starts at day 1 (not zero) --
        dateInds = np.arange(1, len(perfThisSubject)+1)
        # -- Mask out missing sessions (which appear as Nan) --
        notNaN = ~np.isnan(perfThisSubject)  # Used to mask NaN
        # -- Linear fit --
        slope, intercept, rval, pval, se = stats.linregress(dateInds[notNaN],
                                                            perfThisSubject[notNaN])
        slopes[inds] = slope
        intercepts[inds] = intercept

    dframeFit = pd.DataFrame(index=dframe.index)
    dframeFit['slope'] = slopes
    dframeFit['intercept'] = intercepts
    dayToTest = 21
    dframeFit[f'perfAt{dayToTest}d'] = (dayToTest)*slopes + intercepts
    xPercent = 70
    dframeFit[f'daysTo{xPercent}percent'] = (xPercent/100-intercepts)/slopes
    return dframeFit
        
def fit_learning_curves_old(dframePerf, sessionsRange):
    """
    Find linear fit of learning curves.

    Args:
        dframePerf (pd.DataFrame): pandas dataframe containing fraction correct for each session.
        sessionsRange (list): two strings indicating first and last date (in format YYYY-MM-DD)
    Returns:
        dframeSubset (pd.DataFrame): like dframePerf, but excluding bad sessions and subjects
             that did not pass Stage2.
        dframeFit (pd.DataFrame): pandas dataframe with columns related to linear fit
            'slope', 'intercept', 'perfAt21d', 'daysTo70percent'.                       
    """
    #sessionsRange = studyparams.SESSIONS_COH2_STAGE3
    sessionsToPlot = behavioranalysis.sessions_in_range(sessionsRange)
    dframe = dframePerf.copy()
    dframe = dframe[sessionsToPlot]
    
    # -- Exclude mice that took a while not pass criteria 2->3 --
    #miceToExclude = ['pamo020'] #[]#['pamo020','pamo024']
    #dframe = dframe.drop(index=miceToExclude)

    # -- Exclude sessions with likely hardware issues --
    for subject, sessions in studyparams.UNRELIABLE_SESSIONS.items():
        for sessionISO in sessions:
            session = sessionISO.replace('-','')+'a'
            if session in dframe:
                dframe[session][subject] = np.nan
    
    nSubjects = len(dframe)
    slopes = np.empty(nSubjects)
    intercepts = np.empty(nSubjects)

    for inds in range(nSubjects):
        perfThisSubject = dframe.iloc[inds].to_numpy()
        # -- Assumes contiguous dates and starts at day 1 (not zero) --
        dateInds = np.arange(1, len(perfThisSubject)+1)
        # -- Mask out missing sessions (which appear as Nan) --
        notNaN = ~np.isnan(perfThisSubject)  # Used to mask NaN
        # -- Linear fit --
        slope, intercept, rval, pval, se = stats.linregress(dateInds[notNaN],
                                                            perfThisSubject[notNaN])
        slopes[inds] = slope
        intercepts[inds] = intercept

    dframeFit = pd.DataFrame(index=dframe.index)
    dframeFit['slope'] = slopes
    dframeFit['intercept'] = intercepts
    dayToTest = 21
    dframeFit[f'perfAt{dayToTest}d'] = (dayToTest)*slopes + intercepts
    xPercent = 70
    dframeFit[f'daysTo{xPercent}percent'] = (xPercent/100-intercepts)/slopes
    return (dframe, dframeFit)

def find_fast_learners(threshold=21, dframeFit=None):
    """
    Find the group of subjects that reached 70% performance by X days.

    Args:
        threshold (int): number of days to reach 70%.
        dframeFit (pd.DataFrame): optional dataframe with linear fits for stage 3 sessions
    """
    if dframeFit is None:
        # -- Load stage 3 data to find good mice --
        #FIGNAME = 'learning_curve_stage3'
        #figDataFile = 'fraction_correct_stage3.csv'
        #figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
        #figDataFullPath = os.path.join(figDataDir, figDataFile)
        #dframe = pd.read_csv(figDataFullPath, index_col=0)
        #(dframePerf, dframeFit) = fit_learning_curves(dframe)
        dframe = load_stage3(excludeAntibias=1)  # Loads 26 days by default
        dframeFit = fit_learning_curves(dframe)
    fastLearnersBool = dframeFit.daysTo70percent < threshold  # Less than X days to reach 70%
    fastSubjects = list(dframeFit[fastLearnersBool].index)   
    slowSubjects = list(dframeFit[~fastLearnersBool].index)   
    return fastSubjects, slowSubjects
 
def get_slope(curveParams, xsep=0.1):
    """
    Estimate max slope of the psychometric curve.
    """
    slopePointSep = xsep
    slopePointsX = curveParams[0] + slopePointSep*np.array([-0.5, 0.5])
    slopePointsY = extrastats.psychfun(slopePointsX, *curveParams)
    slope = np.diff(slopePointsY)[0]/slopePointSep
    return slope

def plot_psychometric(possibleValues, fractionHitsEachValue, ciHitsEachValue):
    upperWhisker = ciHitsEachValue[1, :]-fractionHitsEachValue
    lowerWhisker = fractionHitsEachValue-ciHitsEachValue[0, :]
    (pline, pcaps, pbars) = plt.errorbar(possibleValues, 100*fractionHitsEachValue,
                                         yerr=[100*lowerWhisker, 100*upperWhisker], color='k')
    pdots, = plt.plot(possibleValues, 100*fractionHitsEachValue, 'o', mec='none', mfc='k', ms=8)
    plt.setp(pline, lw=2)
    ax = plt.gca()
    plt.ylim([0, 100])
    valRange = possibleValues[-1]-possibleValues[0]
    plt.xlim([possibleValues[0]-0.1*valRange, possibleValues[-1]+0.1*valRange])
    return pline, pcaps, pbars, pdots

def session_to_isodate(datestr):
    return datestr[0:4]+'-'+datestr[4:6]+'-'+datestr[6:8]

def days_around(isodate, outformat='%Y-%m-%d'):
    dateformat = '%Y-%m-%d'
    thisDay = datetime.datetime.strptime(isodate, dateformat)
    dayBefore =  thisDay - datetime.timedelta(days=1)
    dayAfter =  thisDay + datetime.timedelta(days=1)
    return (dayBefore.strftime(outformat), thisDay.strftime(outformat),
            dayAfter.strftime(outformat))

def gaussian(x, a, x0, sigma, y0):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+y0

def gaussians_mid_cross(mu1, mu2, std1, std2, amp1, amp2):
    '''
    Return the cross point between the means of the Gaussians
    https://stackoverflow.com/questions/22579434/
    python-finding-the-intersection-point-of-two-gaussian-curves
    '''
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = mu2/(std2**2) - mu1/(std1**2)
    c = mu1**2 /(2*std1**2) - mu2**2 / (2*std2**2) - np.log((std2*amp1)/(std1*amp2))
    allCrosses = np.roots([a,b,c])
    crossInd = np.logical_xor(allCrosses>mu1, allCrosses>mu2)
    midCross = float(allCrosses[crossInd])
    return midCross

def mice_each_condition(group='all', aoc=None):
    """
    group: 'all', 'fast', 'slow'

    aoc (activeOnlyCohorts): 
    activeOnlyCohorts = [2] # Default
    activeOnlyCohorts = [2,3,4]
    """
    # -- Load distributions --
    FIGNAME = 'learning_comparison' 
    figDataDir = os.path.join(settings.FIGURES_DATA_PATH, studyparams.STUDY_NAME, FIGNAME)
    figDataFile = 'distributions_learning.npz'
    figDataFullPath = os.path.join(figDataDir, figDataFile)
    distrib = np.load(figDataFullPath, allow_pickle=True)
    eachCond = list(distrib['eachCond'])
    goodLearners = distrib['goodLearnersT70']
    #miceEachCond = []
    miceEachCond = list(distrib['subjectsEachCond'])
    for indcond, (miceThisCond, goodThisCond) in enumerate(zip(miceEachCond,goodLearners)):
        if group=='fast':
            miceEachCond[indcond] = [mm for mm, good in zip(miceThisCond,goodThisCond) if good]
        elif group=='slow':
            miceEachCond[indcond] = [mm for mm, good in zip(miceThisCond,goodThisCond) if not good]
    # -- Select active mice by cohort --
    if aoc is not None:
        miceFromCohortsToInclude = []
        for indcohort, thisCohort in enumerate(aoc):
            miceFromCohortsToInclude.extend(getattr(studyparams, f'MICE_ALL_COH{thisCohort}'))
            miceEachCond[0] = list(set(miceEachCond[0]).intersection(miceFromCohortsToInclude))
    return (miceEachCond, eachCond)
    

if __name__ == '__main__':
    subject = 'pamo009'; stage = 4
    sessionsList = get_sessions(subject, stage)
    print(sessionsList)
