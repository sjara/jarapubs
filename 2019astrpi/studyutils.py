"""
Additional functions to analyze neural responses to tones and AM sounds.

Based on database_generation_funcs.py (by Devin and Matt)
"""

import numpy as np
import scipy.signal
from scipy import optimize

def cellinfo(celldb, cellID, prefix='', postfix=''):
    """
    cellID can be int or list.
    """
    if isinstance(cellID, int):
        cellID = [cellID]
    for oneCell in cellID:
        if oneCell==-1:
            print('')
        else:
            comment = f'  # cell {oneCell}'
            cc = celldb.loc[oneCell]
            cellStr = f"['{cc.subject}', '{cc.date}', {cc.pdepth}, {cc.egroup}, {cc.cluster}]"
            print(prefix + cellStr + postfix + comment)
        
def select_cells(celldb, restrictND1=False):
    """
    Args:
        restrictND1 (bool): If True, use only nD1 cells that come from tetrodes with D1 cells.
    """
    
    N_FREQ = 16 # HARDCODED
    N_RATE = 11 # HARDCODED

    # -- Exclude recordings from cortex --
    excludedByArea = np.full(len(celldb), False, dtype='bool')
    if 1:
        # celldb.recordingSiteName.unique()
        areasToExclude = ['Supplemental somatosensory area, layer 6b',
                          'Supplemental somatosensory area, layer 6a',
                          'Visceral area, layer 6a',
                          'Primary somatosensory area, barrel field, layer 6b',
                          'Primary somatosensory area, barrel field, layer 6a',
                          'Primary somatosensory area, nose, layer 5']
        for inda, oneArea in enumerate(areasToExclude):
            celldb.recordingSiteName==oneArea
            excludedByArea[celldb.recordingSiteName==oneArea] = True

    # -- Select only responsive cells --
    toneResponsive = celldb['toneMinPval'] < (0.05/N_FREQ)
    amOnsetResponsive = celldb['amMinPvalOnset'] < (0.05/N_RATE)
    amSustainResponsive = celldb['amMinPvalSustain'] < (0.05/N_RATE)
    
    laserBase = celldb.laserBaseRate200
    laserResp = celldb.laserRespRate50
    
    # -- Exclude laser sessions with very low firing rate --
    if 1:
        laserThFR = 1.0 #0.5 # Threshold in spk/s
        highFRlaserSession = (laserBase > laserThFR) | (laserResp > laserThFR)
        #highFRlaserSession = (laserBase > laserThFR)
        print(f'\nWARNING: Excluding cell with laser-session firing rate < {laserThFR} spk/s')
    else:
        highFRlaserSession = True

    cellsWithTone = ~np.isnan(celldb['toneMinPval'])
    cellsWithAM = ~np.isnan(celldb['amMinPvalOnset'])
    cellsWithSoundSessions = cellsWithTone | cellsWithAM
    indD1 = ( (celldb.laserPvalB200R50 < 0.01) & (laserResp>laserBase) &
              highFRlaserSession & ~excludedByArea )
    indND1 = ( (celldb.laserPvalB200R50 > 0.1) | (laserResp<=laserBase) &
               highFRlaserSession & ~excludedByArea )
    #indD1 = (celldb.laserPvalB200R50 < 0.01) & (laserResp>laserBase)
    #indND1 = (celldb.laserPvalB200R50 > 0.05) | (laserResp<=laserBase)
    #indD1 = (celldb.laserPvalB200R50 < 0.05) & (laserResp>laserBase)
    #indND1 = (celldb.laserPvalB200R50 > 0.05) | (laserResp<laserBase)

    # -- Find ND1 in the same tetrode as D1 --
    if restrictND1:
        print('********* WARNING: using only ND1 neighboring D1 ***********')
        celldb['neighborToD1'] = False
        celldbD1 = celldb[indD1]
        otherTetrodeMap = {1:2, 2:1, 3:4, 4:3, 5:6, 6:5, 7:8, 8:7}
        for indRow, dbRow in celldb[indND1].iterrows():
            cells = celldbD1.query('subject==@dbRow.subject and date==@dbRow.date and ' + 
                                   'pdepth==@dbRow.pdepth and egroup==@dbRow.egroup')
            #otherTetrode = otherTetrodeMap[dbRow.egroup]
            #cells = celldbD1.query('subject==@dbRow.subject and date==@dbRow.date and ' + 
            #                       'pdepth==@dbRow.pdepth and ' +
            #                       '(egroup==@dbRow.egroup or egroup==@otherTetrode)')
            if len(cells)>0:
                celldb.at[indRow, 'neighborToD1'] = True
        #print(np.sum(celldb.neighbor & indND1))
        indND1 = indND1 & celldb.neighborToD1
    # -- Find D1 in the same tetrode as ND1 --
    if 0: #restrictD1
        print('********* WARNING: using only D1 neighboring ND1 ***********')
        celldb['neighborToND1'] = False
        celldbND1 = celldb[indND1]
        for indRow, dbRow in celldb[indD1].iterrows():
            cells = celldbND1.query('subject==@dbRow.subject and date==@dbRow.date and ' + 
                                    'pdepth==@dbRow.pdepth and egroup==@dbRow.egroup')
            if len(cells)>0:
                celldb.at[indRow, 'neighborToND1'] = True
        #print(np.sum(celldb.neighbor & indND1))
        indD1 = indD1 & celldb.neighborToND1
        indND1 = indND1 & celldb.neighborToD1

    # -- Select only high enough firing rate --
    if 1:
        thFR = 1.0# 1.0 # Threshold in spk/s
        highFR = ( (celldb.toneFiringRateBaseline > thFR) | (celldb.toneFiringRateBest > thFR) |
                   (celldb.amFiringRateBaseline > thFR) | (celldb.amFiringRateBestOnset > thFR) |
                   (celldb.amFiringRateBestSustain > thFR))
        indD1 = indD1 & highFR
        indND1 = indND1 & highFR
        print(f'\nWARNING: From the total, we analyze only cells with firing rate > {thFR} spk/s\n')
    
    return (toneResponsive, amOnsetResponsive, amSustainResponsive, indD1, indND1)

    

    
def calculate_fra(avgFiringRateRespMap, avgFiringRateBaseline, responseFraction=0.2):
    """
    Calculate the frequency response area.
    The response threshold is defined to be the baseline firing rate plus some
    fraction (responseFraction) of the difference between the baseline and the cell's
    maximum firing rate under any condition.
    (Sutter and Schreiner, https://doi.org/10.1152/jn.1991.65.5.1207)

    Args:
        avgFiringRateRespMap (numpy.ndarray): Average response for each frequency-intensity
            combination. The array should be size (nIntensities, nFrequencies).
        avgFiringRateBaseline (float): Average spontaneous firing rate.
        fraThreshold (float): Firing rate threshold to estimate the frequency response area.

    Returns:
        fra (numpy.ndarray): frequency response area. Boolean array of shape (nInten, nFreq).
        responseThreshold (float): threshold used to define FRA.
    """
    responseThreshold = avgFiringRateBaseline + responseFraction * (avgFiringRateRespMap.max() -
                                                                    avgFiringRateBaseline)
    fra = avgFiringRateRespMap > responseThreshold
    return fra, responseThreshold


def calculate_intensity_threshold(fra, thresholdCF=0.85):
    """, avgFiringRateRespMap
    Calculates intensity threshold index and characteristic frequency (CF) index.
    Neuron's CF is defined as the frequency with the lowest sound intensity inside
    the FRA where 85% (theshold) of the intensities above were also within the FRA.
    
    Args:
        fra (numpy.ndarray): boolean array of shape (nInten, nFreq).
        avgFiringRateRespMap (numpy.ndarray): Average response for each frequency-intensity
            combination. The array should be size (nIntensities, nFrequencies).
        thresholdCF (float): At least this proportion of the intensities above
            must have a response for the freq to be cf.

    Returns:
        intensityInd (int): intensity threshold index
        freqInd (int): characteristic frequency (CF) index
    """
    rows, cols = fra.nonzero()
    lowestIntensityInd = None
    charactFreqInd = None
    cfPoints = []
    for row, col in zip(rows, cols):
        colAbove = fra[row:, col]
        if colAbove.mean() > thresholdCF:
            lowestIntensityInd = row
            charactFreqInd = col
            cfPoints.append((lowestIntensityInd, charactFreqInd))
    if len(cfPoints) == 0:
        lowestIntensityInd, charactFreqInd = None, None
    else:
        cfPointsArray = np.array(cfPoints)
        minPoint = np.argmin(cfPointsArray[:,0])
        lowestIntensityInd, charactFreqInd = cfPointsArray[minPoint,:]
    return lowestIntensityInd, charactFreqInd


def OLD_calculate_intensity_threshold(fra, thresholdCF=0.85):
    """
    Note that there could be more than one frequency that satisfies that criterion;
    this function returns the first one it finds (usually a lower frequency).
    """
    rows, cols = fra.nonzero()
    lowestIntensityInd = None
    charactFreqInd = None
    for row, col in zip(rows, cols):
        colAbove = fra[row:, col]
        if colAbove.mean() > thresholdCF:
            lowestIntensityInd = row
            charactFreqInd = col
            break
    return lowestIntensityInd, charactFreqInd


def calculate_fra_slopes(fra, possibleIntensity, possibleFreq, cfInd, thresholdInd,
                         maxIntensityInd=None):
    """
    Args:
        fra (numpy.ndarray): boolean array of shape (nInten, nFreq).
        possibleIntensity
        possibleFreq
        cfInd (int): characteristic frequency index.
        thresholdInd (int): intensity threshold index.
        maxIntensityInd (int): index of the top of the FRA. Usually calculated as:
                               np.argmax(firingRateRespMap.sum(axis=1))
    Returns:
        lowSlope (float): slope (in dB/oct) of the low-freq side of the FRA
        highSlope (float): slope (in dB/oct) of the high-freq side of the FRA
    """
    ifra = fra.astype(int)
    if maxIntensityInd is None:
        maxIntInd = np.argmax(possibleIntensity)
    else:
        maxIntInd = maxIntensityInd
    medianFiltered = scipy.signal.medfilt(ifra[maxIntInd,:], 3)
    #print(medianFiltered)
    nonZeroInds = np.flatnonzero(medianFiltered)
    if len(nonZeroInds) == 0:
        return (np.nan, np.nan)
    lowBorderInd = nonZeroInds[0]
    highBorderInd = nonZeroInds[-1]
    possibleFreqInOct = np.log2(possibleFreq)
    if lowBorderInd != cfInd:
        lowSlope = ((possibleIntensity[maxIntInd] - possibleIntensity[thresholdInd]) /
                    (possibleFreqInOct[lowBorderInd] - possibleFreqInOct[cfInd]))
    else:
        lowSlope = -np.inf
    if highBorderInd != cfInd:
        highSlope = ((possibleIntensity[maxIntInd] - possibleIntensity[thresholdInd]) /
                     (possibleFreqInOct[highBorderInd] - possibleFreqInOct[cfInd]))
    else:
        highSlope = np.inf
    return (lowSlope, highSlope)
    
def gaussian(x, a, x0, sigma, y0):
    """
    Gaussian function
    Args:
        x (numpy.ndarray): input data
        a (float): the height of the curve's peak
        x0 (float): horizontal offset
        sigma (float): standard deviation (width of Gaussian)
        y0 (float): vertical offset

    Returns:
        output of gaussian function
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2))+y0


def inverse_gaussian(y, a, x0, sigma, y0):
    """
    Inverse function of Gaussian (not to confuse with inverse gaussian distribution).
    This function is useful for finding the lower and upper frequencies
    at the response threshold, given a Gaussian fit at a particular intensity.

    Args:
        y (float or np.array): output value of the Gaussian
        a (float): the height of the curve's peak
        x0 (float): horizontal offset
        sigma (float): standard deviation (width of Gaussian)
        y0 (float): vertical offset

    Returns:
        lower (float): Value of the lower boundary of the Gaussian
        upper (float): Value of the upper boundary of the Gaussian
    """
    sqrtInner = -1*np.log((y-y0)/a)*2*sigma**2
    if sqrtInner < 0:  # No solutions
        return None
    else:
        lower = x0 - np.sqrt(sqrtInner)
        upper = x0 + np.sqrt(sqrtInner)
        return [lower, upper]

'''
# -- Testing the inverse Gaussian --
plt.clf()
x = np.arange(-10,10,0.1)
a=2; x0=1; sigma=2; y0=1;
y = a*np.exp(-(x-x0)**2/(2*sigma**2))+y0
plt.plot(x,y,'.');

xsq = -2*(sigma**2) * np.log((y-y0)/a)
x = (x0 + np.sqrt(xsq))
plt.plot(x,y,'.');
'''


def calculate_BWx(intensityEachTrial, freqEachTrial, firingRateEachTrial,
                   intensityThreshold, dBabove):
    # GAUSSIAN PARAMS: a, x0, sigma, y0
    possibleFreq = np.unique(freqEachTrial)
    possibleIntensity = np.unique(intensityEachTrial)
    possibleLogFreq = np.log2(possibleFreq)
    intensityIndXdBabove = np.argmin(np.abs(possibleIntensity-intensityThreshold-dBabove))
    selectedTrialsByIntensity = intensityEachTrial==possibleIntensity[intensityIndXdBabove]
    firingRateSelIntensity = firingRateEachTrial[selectedTrialsByIntensity]
    logFreqsSelIntensity = np.log2(freqEachTrial[selectedTrialsByIntensity])
    p0 = [1, possibleLogFreq[len(possibleLogFreq)//2], 1, np.min(firingRateSelIntensity)]
    bounds = ([0, possibleLogFreq[0], 0, 0],
              [np.inf, possibleLogFreq[-1], np.inf, np.inf])
    try:
        #popt, pcov = optimize.curve_fit(gaussian, logFreqsSelIntensity,
        #                                firingRateSelIntensity, p0=p0, bounds=bounds)
        popt, pcov = optimize.curve_fit(gaussian, logFreqsSelIntensity,
                                        firingRateSelIntensity, p0=p0)
    except RuntimeError:
        print("Could not fit gaussian curve to tuning data.")
        popt = None
        Rsquared = np.nan
        fullWidthHalfMax = np.nan
    else:
        if popt[2]<0:
            print('Warning! Gaussian fit returned negative standard deviation. Changed to positive.')
            popt[2] = -popt[2]
        gaussianResp = gaussian(logFreqsSelIntensity, *popt)
        residuals = firingRateSelIntensity - gaussianResp
        ssquared = np.sum(residuals**2)
        ssTotal = np.sum((firingRateSelIntensity-np.mean(firingRateSelIntensity))**2)
        Rsquared = 1 - (ssquared/ssTotal)
        fullWidthHalfMax = 2.355*popt[2] # Sigma is popt[2]
    #breakpoint()
    return (popt, Rsquared, fullWidthHalfMax)
   

'''
lowerFreq, upperFreq, Rsquared10AboveSIT = funcs.calculate_BW10_params(ind10Above, popts, Rsquareds,
                                                                       responseThreshold,
                                                                       intensityThreshold)
                if (lowerFreq is not None) and (upperFreq is not None):
                    fitMidpoint = np.sqrt(lowerFreq * upperFreq)
                    bw10 = (upperFreq - lowerFreq) / cf
'''




def angle_population_vector_zar(angles):
    """
    Computes the length of the mean vector for a population of angles.
    Copied directly from Biostatistical analysis, Zar, 3rd ed, pg 598
    (Mike Wehr has this book)
    Args:
        angles (numpy.ndarray): Each value is an angle in radians
    Returns:
        r (float): Angle population vector for calculations in Rayleigh
        test
    """
    X = np.mean(np.cos(angles))
    Y = np.mean(np.sin(angles))
    r = np.sqrt(X**2 + Y**2)
    return r


def rayleigh_test(angles):
    """
    Performs Rayleigh Test for non-uniformity of circular data.
    Compares against Null hypothesis of uniform distribution around circle
    Assume one mode and data sampled from Von Mises.
    Use other tests for different assumptions.
    Maths from [Biostatistical Analysis, Zar].
    Args:
        angles (numpy.ndarray): Each value is an angle in radians
    Returns:
        zVal (float): Statistic from Rayleigh test
        pVal (float): Significance value from Rayleigh test
    """
    N = len(angles)
    # Compute Rayleigh's R
    R = N * angle_population_vector_zar(angles)
    # Compute Rayleight's z
    zVal = (R**2) / N
    # Compute pvalue (Zar, Eq 27.4)
    pVal = np.exp(np.sqrt(1. + 4*N + 4*(N**2. - R**2)) - 1. - 2.*N)
    return zVal, pVal

