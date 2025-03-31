"""
This file contains the default names and queries used for this study.
"""

from jaratoolbox import colorpalette as cp
import matplotlib


STUDY_NAME = '2023acid'

SUBJECTS = ['acid006', 'acid007', 'acid009', 'acid010']
REAGENTS = ['pre', 'saline', 'doi']

# Date when pre/post sync light was added to videos.
PREPOST_SYNC_LIGHT_START_DATE = '2023-08-08'

# Date of session with no sync light.
N0_LIGHT_DATE = '2023-03-22'

N_FREQ = 16

STIMULI = ['high', 'low', 'down', 'up']

ODDBALL_SESSION_TYPES = {'HighFreq': {'video':'01', 'behavior':'a', 'paradigm':'oddball_sequence',
                                      'oddball':'High', 'standard':'Low'},
                         'LowFreq': {'video':'02', 'behavior':'b', 'paradigm':'oddball_sequence',
                                     'oddball':'Low', 'standard':'High'},
                         'FM_Down': {'video':'03', 'behavior':'c', 'paradigm':'oddball_sequence',
                                     'oddball':'Down', 'standard':'Up'},
                         'FM_Up': {'video':'04', 'behavior':'d', 'paradigm':'oddball_sequence',
                                   'oddball':'Up', 'standard':'Down'},
                         }

RUNNING_THRESHOLD = 3

MIN_R_SQUARED = 0.05 #0.1 #0.05

MIN_SIGMA = 0.15  # Minimum value of Sigma for the Gaussian fit.
MAX_SIGMA = 6     # Maximum value of Sigma for the Gaussian fit.

FR_THRESHOLD = 5  # Minimum evoked firing rate to consider a cell responsive.

MAX_CHANGE_FACTOR = 1.3  # Maximum change factor to consider a cell steady

