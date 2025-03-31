"""
Plot examples of spike templates.
"""

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from jaratoolbox import extraplots
from jaratoolbox import celldatabase
from jaratoolbox import settings
import studyparams
import figparams
import importlib
importlib.reload(figparams)

'''
# -- Find the cells we'll plot --
dbPath = os.path.join(settings.DATABASE_PATH, studyparams.STUDY_NAME)
dbFilename = os.path.join(dbPath, f'celldb_{studyparams.STUDY_NAME}_basic.h5')
celldb = celldatabase.load_hdf(dbFilename)
celldb[(celldb.subject=='acid010') & (celldb.date=='2023-08-07') & (celldb.cluster==34)]  # 1133
celldb[(celldb.subject=='acid010') & (celldb.date=='2023-08-07') & (celldb.cluster==36)]  # 1135
'''

subject = 'acid010'
dataDir = settings.EPHYS_NEUROPIX_PATH
multiDir = 'multisession_2023-08-07_3000um_processed'
clusterFolder = os.path.join(dataDir, subject, multiDir)

clusterToPlot = [34,36]
colorEachCluster = ['teal', 'olivedrab']
centerChan = 152;
ylimOffset = 80

tscale = 0.25
ascale = 100

clusterGroup = pd.read_csv(os.path.join(clusterFolder, 'cluster_group.tsv'), sep='\t')
amplitude = np.load(os.path.join(clusterFolder, 'amplitudes.npy')).squeeze()
clusters = np.load(os.path.join(clusterFolder, 'spike_clusters.npy')).squeeze()
templates = np.load(os.path.join(clusterFolder, 'templates.npy'))
(nOrigClusters, nTimePoints, nChannels) = templates.shape
channelPos = np.load(os.path.join(clusterFolder, 'channel_positions.npy'))
channelMap = np.load(os.path.join(clusterFolder, 'channel_map.npy'))

goodClusters = clusterGroup.cluster_id[clusterGroup['group']=='good'].to_numpy()

def plot_template(oneTemplate, channelPos, ascale=100, tscale=0.25,
                  subset=True, centerChan=None, color='k'):
    if subset:
        indMax = np.argmax(np.abs(oneTemplate))
        (sampleMax, chanMax) = np.unravel_index(indMax, oneTemplate.shape)
        if centerChan is not None:
            chanMax = centerChan
        channelsToPlot = np.arange(chanMax-6-chanMax%2, chanMax+8-chanMax%2)
    else:
        channelsToPlot = np.arange(nChannels)
    tvec = tscale * np.arange(nTimePoints)
    samplesToPlot = np.arange(21, nTimePoints)
    for oneChannel in channelsToPlot:
        xoffset = channelPos[oneChannel,0]
        yoffset = channelPos[oneChannel,1]
        plt.plot(tvec[samplesToPlot]+xoffset,
                 ascale*oneTemplate[samplesToPlot,oneChannel]+yoffset, color=color)

plt.clf()
for indc, oneCluster in enumerate(clusterToPlot):
    indsThisCluster = (clusters==oneCluster)
    medianAmp = np.median(amplitude[indsThisCluster])
    oneTemplate = templates[oneCluster,:,:]
    indMax = np.argmax(np.abs(oneTemplate))
    (sampleMax, chanMax) = np.unravel_index(indMax, oneTemplate.shape)
    plt.subplot(1,2,1+indc)
    plot_template(oneTemplate, channelPos, centerChan=centerChan, color=colorEachCluster[indc])
    plt.title(f'Cluster {oneCluster}, ch={chanMax}\n(median amplitude {medianAmp:0.2f})')
    plt.ylim(channelPos[centerChan,1]-ylimOffset, channelPos[centerChan,1]+ylimOffset)
plt.show()

extraplots.save_figure('plots_spike_templates', 'svg', [4,4], outputDir='/tmp/')


