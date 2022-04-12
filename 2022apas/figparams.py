"""
Common parameters for figures related to this study.
"""

from jaratoolbox import colorpalette as cp
import matplotlib

matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none'  # To render as font rather than outlines

STUDY_NAME = '2022apas'

fontSizeLabels = 12
fontSizeTicks = 12
fontSizePanel = 16

colors = {}
colors['activeOnly'] = cp.TangoPalette['SkyBlue2']
colors['activePassive'] = cp.TangoPalette['Chameleon3']

