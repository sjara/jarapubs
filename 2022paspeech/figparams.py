"""
Common parameters for figures related to this study.
"""

from jaratoolbox import colorpalette as cp
import matplotlib

#matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none'  # To render as font rather than outlines

STUDY_NAME = '2022paspeech'

fontSizeLabels = 16
fontSizeTicks = 12
fontSizePanel = 20

colors = {}
#colors['blueLaser'] = cp.TangoPalette['SkyBlue1']
colors['sound'] = cp.TangoPalette['Butter2']
colors['ft'] = cp.TangoPalette['ScarletRed1']
colors['vot'] = cp.TangoPalette['SkyBlue2']
colors['audP'] = cp.TangoPalette['SkyBlue1']
colors['audD'] = cp.TangoPalette['Orange1']
colors['audV'] = cp.TangoPalette['Aluminium4']
