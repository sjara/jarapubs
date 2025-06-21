"""
Additional functions for 2024traje project.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

def simplify_site_name(recordingSiteName):
    simplerName = recordingSiteName.split(',')[0]
    simplerName = simplerName.replace('Supplemental', 'Supp')
    return simplerName

def average_spike_count(spikeCount, trialsEachCond):
    """
    Args:
        spikeCounts (np.ndarray): [nInstances, nTrials, nTimebins]
        trialsEachCond (np.ndarray): [nTrials, nConditions]
    Returns:
        avgSpikeCountsEachCond (np.ndarray): [nCells, nInstances, nTimeBins]
    """
    nCells, nTrials, nTimebins = spikeCount.shape
    nInstances = trialsEachCond.shape[1]
    avgSpikeCountsEachCond = np.empty([nCells, nInstances, nTimebins])
    for indcond in range(nInstances):
        spikeCountsThisCond = spikeCount[:, trialsEachCond[:,indcond], :]
        avgSpikeCountsEachCond[:, indcond, :] = spikeCountsThisCond.mean(axis=1)
    return avgSpikeCountsEachCond

def smooth_spike_count(spikeCounts, winsize=5, winfun=np.hanning):   
    """
    Args:
        spikeCounts (np.ndarray): [nCells, nConditions, nTimebins]
        winsize (int): size of smoothing window
        winfun (function): smoothing window
    Returns:
        smoothSpikeCounts (np.ndarray): [nCells, nConditions, nTimebins]
    """
    window = winfun(winsize)
    window /= window.sum()
    smoothSpikeCounts = np.empty_like(spikeCounts)
    smoothingFunction = lambda m: np.convolve(m, window, mode='same')
    for indcond in range(smoothSpikeCounts.shape[1]):
        smoothSpikeCounts[:, indcond, :] = \
           np.apply_along_axis(smoothingFunction, axis=1, arr=spikeCounts[:, indcond, :]) 
    return smoothSpikeCounts

def euclidian_distance(v1, v2):
    """v1 and v2 are 2D arrays of shape (nCells, nTimebins)"""
    return np.linalg.norm(v1-v2, axis=0)

def cosine_distance(v1, v2):
    """v1 and v2 are 2D arrays of shape (nCells, nTimebins)"""
    dot_product = np.sum(v1 * v2, axis=0)
    v1_norm = np.linalg.norm(v1, axis=0)
    v2_norm = np.linalg.norm(v2, axis=0)
    cosine_similarity = dot_product / (v1_norm * v2_norm)
    return 1 - cosine_similarity

def distance_between_trajectories(spikeCount, distanceType):
    """
    Args:
        spikeCounts (np.ndarray): [nCells, nConditions, nTimebins]
        distanceType (str): 'cosine' or 'euclidian'
    Returns:
        distance (np.ndarray): [nConditions, nConditions, nTimebins]
    """
    if distanceType=='euclidian':
        distancefun = euclidian_distance
    elif distanceType=='cosine':
        distancefun = cosine_distance
    else:
        raise ValueError('distanceType must be either euclidian or cosine')
    nCells, nConditions, nTimebins = spikeCount.shape
    distance = np.empty([nConditions, nConditions, nTimebins])
    for indcond1 in range(nConditions):
        for indcond2 in range(nConditions):
            if indcond1==indcond2:
                distance[indcond1, indcond2, :] = 0
                continue
            if indcond1>indcond2:
                distance[indcond1, indcond2, :] = distance[indcond2, indcond1, :]
                continue
            distance[indcond1, indcond2, :] = distancefun(spikeCount[:, indcond1, :],
                                                          spikeCount[:, indcond2, :])
    return distance


def plot_distance_matrix(distance, binEdges, periods, vlims=[None,None], fig=None):
    """
    Args:
        distance (np.ndarray): [nConditions, nConditions, nTimebins]
        periods (list): timebins to plot
    """
    nPeriods = len(periods)
    if fig is None:
        fig = plt.gcf()
    axs = fig.subplots(2, nPeriods//2, sharex=True, sharey=True)
    for indp in range(nPeriods):
        sampleThisPeriod = np.argmin(np.abs(binEdges-periods[indp]))
        im = axs.flat[indp].imshow(distance[:,:,sampleThisPeriod], aspect='equal',
                                   vmin=vlims[0], vmax=vlims[1])
        axs.flat[indp].set_title(f't = {periods[indp]}')

        
def get_onset_offset_bins(spkData):
    """
    Args:
        spkData (dict): must contain keys 'binEdges', 'soundDuration'
    """
    onsetBin = np.flatnonzero(spkData['binEdges']<0)[-1]
    offsetBin = np.flatnonzero(spkData['binEdges']<spkData['soundDuration'])[-1]
    return onsetBin, offsetBin

def plot_distance_traces(distance, spkData, refInstance=0, xlim=None, ylim=None):
    """
    Args:
        distance (np.ndarray): [nConditions, nConditions, nTimebins]
        spkData (dict): must contain keys 'trialsEachInstance', 'trialsEachCateg', 'binEdges'
        refInstance (int): reference instance to compare to all others
    """
    nCategories = spkData['trialsEachCateg'].shape[1]
    nInstances = spkData['trialsEachInstance'].shape[1]
    nInstancesPerCategory = nInstances//nCategories
    timeBins = spkData['binEdges'][:-1]
    fig = plt.gcf()
    gs = gridspec.GridSpec(nCategories, nCategories, figure=fig, wspace=0.05, hspace=0.2)
    axs = []
    #axs = [fig.add_subplot(gs[ind]) for ind in range(nCategories)]
    #axs = [fig.add_subplot(gs[ind], sharex=None if ind == 0 else axs[0],
    #                       sharey=None if ind == 0 else axs[0]) for ind in range(nCategories)]
    colorEachCategory = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for indref in range(nCategories):
        for indcomp in range(nCategories):
            ax = fig.add_subplot(gs[indref, indcomp])
            ax.axvline(0, color='0.5', linestyle='-')
            refInstance = indref*nInstancesPerCategory
            compInstances = np.arange(indcomp*nInstancesPerCategory, (indcomp+1)*nInstancesPerCategory)
            thisColor = colorEachCategory[compInstances[0]//nInstancesPerCategory]
            ax.plot(timeBins, distance[refInstance, compInstances, :].T, color=thisColor)
            ax.set_title(f'Instance {refInstance} to {compInstances}')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.grid(True)
            ax.set_ylabel('Distance')
            ax.set_xlabel('Time (s)')
            if indcomp==indref:
                ax.set_facecolor('0.95')
            ax.label_outer()
            axs.append(ax)
    return axs


def plot_trajectories(avgSpikeCountsEachCond, onsetBin, offsetBin, labels, comps=[0,1], verb=False):
    """
    Args:
        avgSpikeCountsEachCond (np.ndarray): [nCells, nInstances, nTimeBins]
        onsetBin (int): index of sound onset timebin
        offsetBin (int): index of sound offset timebin
        labels (list): list of labels for each condition
        comps (list): list of principal components (or cells) to plot, starting from 0
    """
    nCategories = avgSpikeCountsEachCond.shape[1]
    PCx, PCy = comps
    nRows = 2
    nCols = int(np.ceil(nCategories/nRows))
    gsMain = gridspec.GridSpec(nRows, nCols)
    gsMain.update(left=0.05, right=0.98, top=0.95, bottom=0.06, wspace=0.25, hspace=0.3)
    axs = [plt.subplot(gsMain[i, j]) for i in range(nRows) for j in range(nCols)]
    axAll = axs[-1] #plt.subplot(gsBottom[nRows-1, nCols-1])
    colorEachCategory = plt.rcParams['axes.prop_cycle'].by_key()['color'] # blue, orange, green, red, purple
    alphaColor = 0.5
    colorBaseline = 'gray'
    markerSize = 4
    lineWidth = 1
    for indcond in range(nCategories):
        ax = axs[indcond]
        if indcond>0:
            ax.sharex(axs[0])
            ax.sharey(axs[0])
        ax.plot(avgSpikeCountsEachCond[PCx, indcond, :onsetBin+1],
                avgSpikeCountsEachCond[PCy, indcond, :onsetBin+1],
                '.-', color=colorBaseline, alpha=alphaColor, ms=markerSize, lw=lineWidth)
        ax.plot(avgSpikeCountsEachCond[PCx, indcond, onsetBin:offsetBin+1],
                avgSpikeCountsEachCond[PCy, indcond, onsetBin:offsetBin+1],
                '.-', color=colorEachCategory[indcond], alpha=alphaColor, ms=markerSize, lw=lineWidth)
        ax.plot(avgSpikeCountsEachCond[PCx, indcond, offsetBin:],
                avgSpikeCountsEachCond[PCy, indcond, offsetBin:],
                '.-', color=colorBaseline, alpha=alphaColor, ms=markerSize, lw=lineWidth)
        ax.plot(avgSpikeCountsEachCond[PCx, indcond, onsetBin],
                avgSpikeCountsEachCond[PCy, indcond, onsetBin],
                'o', ms=4, color='g')
        ax.plot(avgSpikeCountsEachCond[PCx, indcond, offsetBin],
                avgSpikeCountsEachCond[PCy, indcond, offsetBin],
                'o', ms=4, color='r')
        if indcond//3==1:
            ax.set_xlabel(f'PC {PCx+1}')
        if indcond%3==0:
            ax.set_ylabel(f'PC {PCy+1}')
        ax.set_title(labels[indcond], fontweight='bold')

        axAll.plot(avgSpikeCountsEachCond[PCx, indcond], avgSpikeCountsEachCond[PCy, indcond],
                    '-', color=colorEachCategory[indcond], alpha=0.5)
        axAll.sharex(axs[0])
        axAll.sharey(axs[0])
    if verb:
        print('LINES: Gray: baseline. Color: evoked')
        print('DOTS: Green: onset. Red: offset.')
                      

def old_plot_distance_traces(distance, toCompare, spkData, ncols=2):
    """
    Args:
        distance (np.ndarray): [nConditions, nConditions, nTimebins]
        toCompare (list): list of lists, each containing a reference instance and a list of instances to compare
        spkData (dict): must contain keys 'trialsEachInstance', 'trialsEachCateg', 'binEdges'
    """
    nInstancesPerCategory = spkData['trialsEachInstance'].shape[1]//spkData['trialsEachCateg'].shape[1]
    timeBins = spkData['binEdges'][:-1]
    fig = plt.gcf()
    gs = gridspec.GridSpec(len(toCompare)//ncols, ncols, figure=fig, wspace=0.05, hspace=0.2)
    axs = [fig.add_subplot(gs[ind]) for ind in range(len(toCompare))]
    colorEachCategory = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for indcomp, comparison in enumerate(toCompare):
        ax = axs[indcomp]
        refInstance, compInstances = comparison
        thisColor = colorEachCategory[compInstances[0]//nInstancesPerCategory]
        ax.plot(timeBins, distance[refInstance, compInstances, :].T, color=thisColor)
        ax.set_title(f'Instance {refInstance} to {compInstances}')
        ax.grid(True)
        ax.set_ylabel('Distance')
        ax.set_xlabel('Time (s)')
        ax.label_outer()
    return axs




'''
toCompare = [ [0,  [ 0,  1,  2,  3]],
              [0,  [12, 13, 14, 15]],
              [0,  [16, 17, 18, 19]],
              [12, [ 0,  1,  2,  3]],
              [12, [12, 13, 14, 15]],
              [12, [16, 17, 18, 19]],
              [16, [ 0,  1,  2,  3]],
              [16, [12, 13, 14, 15]],
              [16, [16, 17, 18, 19]] ]
'''
'''
    # The next line doesn't work on my version of matplotlib
    #axs = fig.subplots(1, 4, gridspec_kw={'wspace': 0.2}, sharex=True, sharey=True)
'''
