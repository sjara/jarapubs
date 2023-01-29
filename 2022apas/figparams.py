"""
Common parameters for figures related to this study.
"""

from jaratoolbox import colorpalette as cp
import matplotlib

matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none'  # To render as font rather than outlines

fontSizeLabels = 8
fontSizeTicks = 8
fontSizePanel = 12

colors = {}
colors['activeOnly'] = cp.TangoPalette['SkyBlue2']
colors['activePassive'] = cp.TangoPalette['Chameleon3']
#colors['passiveThenActive'] = cp.TangoPalette['Orange3']
colors['passiveThenActive'] = cp.TangoPalette['ScarletRed2']
'''

colors['activeOnly'] = cp.TangoPalette['ScarletRed2']
colors['activePassive'] = cp.TangoPalette['SkyBlue2']
colors['passiveThenActive'] = cp.TangoPalette['Chameleon3']
'''

