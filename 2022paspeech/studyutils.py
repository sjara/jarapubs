"""
Additional functions to analyze neural responses to tones and AM sounds.

Based on database_generation_funcs.py (by Devin and Matt)
"""

import numpy as np

# -- Methods to convert from pixels of Atlas to mm --
def pix2mmAP(pixels):
    return -0.94 - (280-pixels)*0.025
def mm2pixAP(mm):
    return ((mm+0.94) / 0.025) + 280
def pix2mmDV(pixels):
    return (pixels-10)*0.025
def mm2pixDV(mm):
    return (mm/0.025) + 10


LATEX_TABLE_RESPONSIVE_HEADER = """
\\begin{tabular}{|c|c|c|}
\\hline
\\textbf{Area} & \\textbf{Percentage} & \\textbf{Responsive/Total} \\\\ \\hline \
"""
LATEX_TABLE_RESPONSIVE_DATAROW = """
{area} & {percent:0.1f}\\% & {num}/{den} \\\\ \\hline \
"""
LATEX_TABLE_RESPONSIVE_FOOTER = """
\\end{tabular}
"""

LATEX_TABLE_PVALS_HEADER = """
\\begin{tabular}{|c|c|c|c|c|}
\\hline \
"""
LATEX_TABLE_PVALS_DATAROW = """
\\textbf{{{area}}} & {0} & {1} & {2} & {3} \\\\ \\hline \
"""
LATEX_TABLE_PVALS_FOOTER = """
\\end{tabular}
"""

LATEX_TABLE_SELECTIVE_HEADER ="""
\\begin{tabular}{|c|c|c|c|c|}
\\hline
\\textbf{Area} & \\textbf{Of total} & \\textbf{Sel/Total} & \\textbf{\\makecell{Of sound\\\\[-1ex] responsive}} & \\textbf{\\makecell{Of speech\\\\[-1ex] responsive}} \\\\ \\hline \
"""
LATEX_TABLE_SELECTIVE_DATAROW ="""
{area} & {pcOfTotal:0.1f}\\% & {num}/{denomsTotal} & {pcOfSoundResp:0.1f}\\% & {pcOfSpeechResp:0.1f}\\% \\\\ \\hline \
"""
LATEX_TABLE_SELECTIVE_FOOTER ="""
\end{tabular}
"""

LATEX_TABLE_MIXSELECTIVE_HEADER ="""
\\begin{tabular}{|c|c|c|c|c|}
\\hline
\\textbf{Area} & \\textbf{Of total} & \\textbf{Sel/Total} & \\textbf{\\makecell{Of sound\\\\[-1ex] responsive}} & \\textbf{\\makecell{Of speech\\\\[-1ex] responsive}} \\\\ \\hline \
"""
LATEX_TABLE_MIXSELECTIVE_DATAROW ="""
{area} & {pcOfTotal:0.1f}\\% & {num}/{denomsTotal} & {pcOfSoundResp:0.1f}\\% & {pcOfSpeechResp:0.1f}\\% \\\\ \\hline \
"""
LATEX_TABLE_MIXSELECTIVE_FOOTER ="""
\end{tabular}
"""


def latex_table_mixselective(areas, numerators, denomsTotal, denomsSoundResp, denomsSpeechResp,
                             parcellation='Area'):
    """
    Args:
        parcellation (str): Name of parcellation, e.g., 'Area' or 'Region'.
    """
    tableStr = ''
    tableStr += LATEX_TABLE_MIXSELECTIVE_HEADER
    for inda, oneArea in enumerate(areas):
        pcOfTotal = 100 * numerators[inda] / denomsTotal[inda]
        pcOfSoundResp = 100 * numerators[inda] / denomsSoundResp[inda]
        pcOfSpeechResp = 100 * numerators[inda] / denomsSpeechResp[inda]
        tableStr += LATEX_TABLE_MIXSELECTIVE_DATAROW.format(area = oneArea,
                                                         pcOfTotal = pcOfTotal,
                                                         num = numerators[inda],
                                                         denomsTotal = denomsTotal[inda],
                                                         pcOfSoundResp = pcOfSoundResp,
                                                         pcOfSpeechResp = pcOfSpeechResp)
    tableStr += LATEX_TABLE_MIXSELECTIVE_FOOTER
    tableStr.replace('Area', parcellation)
    return tableStr

def latex_table_responsive(areas, numerators, denominators, parcellation='Area'):
    """
    Args:
        parcellation (str): Name of parcellation, e.g., 'Area' or 'Region'.
    """
    tableStr = ''
    tableStr += LATEX_TABLE_RESPONSIVE_HEADER
    for inda, oneArea in enumerate(areas):
        percent = 100 * numerators[inda] / denominators[inda]
        tableStr += LATEX_TABLE_RESPONSIVE_DATAROW.format(area = oneArea,
                                                          percent = percent,
                                                          num = numerators[inda],
                                                          den = denominators[inda])
    tableStr += LATEX_TABLE_RESPONSIVE_FOOTER
    tableStr.replace('Area', parcellation)
    return tableStr


def latex_table_selective(areas, numerators, denomsTotal, denomsSoundResp, denomsSpeechResp,
                          parcellation='Area'):
    """
    Args:
        parcellation (str): Name of parcellation, e.g., 'Area' or 'Region'.
    """
    tableStr = ''
    tableStr += LATEX_TABLE_SELECTIVE_HEADER
    for inda, oneArea in enumerate(areas):
        pcOfTotal = 100 * numerators[inda] / denomsTotal[inda]
        pcOfSoundResp = 100 * numerators[inda] / denomsSoundResp[inda]
        pcOfSpeechResp = 100 * numerators[inda] / denomsSpeechResp[inda]
        tableStr += LATEX_TABLE_SELECTIVE_DATAROW.format(area = oneArea,
                                                         pcOfTotal = pcOfTotal,
                                                         num = numerators[inda],
                                                         denomsTotal = denomsTotal[inda],
                                                         pcOfSoundResp = pcOfSoundResp,
                                                         pcOfSpeechResp = pcOfSpeechResp)
    tableStr += LATEX_TABLE_SELECTIVE_FOOTER
    tableStr.replace('Area', parcellation)
    return tableStr


def pval_to_str(pvalue):
    if np.isnan(pvalue):
        pvalStr = '---'
    elif pvalue < 0.001:
        pvalStr = '$<$0.001'
    else:
        pvalStr = f'{pvalue:0.3f}'.rstrip('0').rstrip('.')
    return pvalStr
    
def latex_table_pvals(areas, pvalues):
    """
    """
    tableStr = ''
    tableStr += LATEX_TABLE_PVALS_HEADER
    titleRow = ' & '.join([f'\\textbf{{{area}}}' for area in ['']+areas]) + ' \\\\ \\hline'
    tableStr += titleRow
    for inda, oneArea in enumerate(areas):
        pvalsList = [pval_to_str(pval) for pval in pvalues[inda]]
        tableStr += LATEX_TABLE_PVALS_DATAROW.format(area = oneArea, *pvalsList)
    tableStr += LATEX_TABLE_PVALS_FOOTER
    return tableStr



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
    """
    Calculates intensity threshold index and characteristic frequency (CF) index.
    Neuron's CF is defined as the frequency with the lowest sound intensity inside
    the FRA where 85% (theshold) of the intensities above were also within the FRA.
    Note that there could be more than one frequency that satisfies that criterion;
    this function returns the first one it finds (usually a lower frequency).

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
    for row, col in zip(rows, cols):
        colAbove = fra[row:, col]
        if colAbove.mean() > thresholdCF:
            lowestIntensityInd = row
            charactFreqInd = col
            break
    return lowestIntensityInd, charactFreqInd


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


def OLD_inverse_gaussian(y, a, x0, sigma, y0):
    """
    Inverse function of Gaussian (not to confuse with inverse gaussian distribution).
    This function finds the lower and upper frequencies
    with desirable response threshold output, and optimal parameters (*popt)
    for 10dB above the current bandwidth

    Args:
        y (numpy.ndarray): Gaussian fit
        a (float): the height of the curve's peak
        x0 (float): the position of the center of the peak
        sigma (float): standard variation that determines the width of the 'bell'
        y0 (float): constant

    Returns:
        lower (flaot): Value of the lower boundary of the gaussian at a set
        intensity
        upper (float): Value of the upper boundary of the gaussian at a set
        intensity
    """
    sqrtInner = -1*np.log((y-y0)/a)*2*sigma**2
    if sqrtInner < 0:  # No solutions
        return None
    else:
        lower = 2**(x0 - np.sqrt(sqrtInner))
        upper = 2**(x0 + np.sqrt(sqrtInner))
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


def scatter_hist(x, y, ax, ax_histx, ax_histy, binwidth = None):
    '''
    Plots scatterplot with marginal histograms. Default binwidth is 0.1
    '''
    if binwidth is None:
        binwidth = 0.1
    xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    lim = (int(xymax/binwidth) + 1) * binwidth
    bins = np.arange(-lim, lim+binwidth, binwidth)
    ax_histx.tick_params(axis = "x", labelbottom = True)
    ax_histy.tick_params(axis = "y", labelleft = True)
    ax.scatter(x, y)
    ax_histx.hist(x, bins = bins)
    ax_histy.hist(y, bins = bins, orientation = 'horizontal')
