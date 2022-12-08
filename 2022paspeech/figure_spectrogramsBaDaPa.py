'''
This figure plots 3 spectrograms: /ba/, /pa/, and /da/ at 8x human freq range
'''
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
from jaratoolbox import soundanalysis
from scipy.io import wavfile
import figparams
import studyparams
from importlib import reload
reload(figparams)

FIGNAME = 'figure_spectrogramsBaDaPa'
SAVE_FIGURE = 1
outputDir = 'C:/Users/jenny/tmp/'
figFilename = 'figure_spectrogramsBaDaPa' # Do not include extension
figFormat = 'pdf' # 'pdf' or 'svg'
figSize = [8.5, 3.5] # In inches

PANELS = [3,1] # Plot panel i if PANELS[i]==1

fontSizeLabels = figparams.fontSizeLabels
fontSizeTicks = figparams.fontSizeTicks
fontSizePanel = figparams.fontSizePanel

labelPosX = [0.07, 0.37, 0.7]   # Horiz position for panel labels
labelPosY = [0.95]    # Vert position for panel labels

plt.figure()
gsMain = gridspec.GridSpec(1, 3)
gsMain.update(left=0.1, right=0.9, top=0.85, bottom=0.15, wspace=0.2, hspace=0.4)

soundsDir = 'H:\\jarasounds\\ft_vot_8x_20220115' #HARDCODED
#soundsDir = 'H:\\jarasounds\\ft_vot_1x_20220115'
votValues = ['000', '100']
ftValues = ['100', '000']

a = 0
for indVOT, thisVOT in enumerate(votValues):
    for indFT, thisFT in enumerate(ftValues):
        if indFT==0 and indVOT:
            pass
        else:
            thisSound = f'syllable_8x_vot{thisVOT}_ft{thisFT}.wav'
            soundFile = os.path.join(soundsDir, thisSound)
            samplingRate, wave = wavfile.read(soundFile)
        #nSamples, samplingRate = soundanalysis.SoundAnalysis.get_sound_parameters(data)
        #timeVec, wave = soundanalysis.SoundAnalysis.get_waveform(data)

            plt.subplot(gsMain[0, a])
            soundanalysis.plot_spectrogram(wave, samplingRate)
            plt.ylim([0,28000])
            #plt.title(thisSound)
            plt.xlim([0,.13])
            plt.xticks(np.arange(0, .13, 0.04),['0', '40', '80', '120'], fontsize=fontSizeTicks)
            plt.xlabel('Time (ms)', fontsize=fontSizeLabels)

            if a == 0:
                plt.ylabel('Frequency (kHz)', fontsize=fontSizeLabels)
                plt.yticks(np.arange(0, 28000, 5000), ['0', '5', '10', '15', '20', '25'], 
                    fontsize=fontSizeTicks)
                plt.title('/da/', fontsize=fontSizeLabels, fontweight='bold',)
            elif a ==1 :
                plt.title('/ba/', fontsize=fontSizeLabels, fontweight='bold')
                plt.yticks([])
                plt.ylabel('')
            elif a == 2:
                plt.title('/pa/', fontsize=fontSizeLabels, fontweight='bold')
                plt.yticks([])
                plt.ylabel('')
            a = a + 1


plt.show()

if SAVE_FIGURE:
    extraplots.save_figure(figFilename, figFormat, figSize, outputDir)
