"""
Common parameters for figures related to this study.
"""

from jaratoolbox import colorpalette as cp
import matplotlib

matplotlib.rcParams['font.family'] = 'Helvetica'
matplotlib.rcParams['svg.fonttype'] = 'none'  # To render as font rather than outlines

STUDY_NAME = '2022paspeech'

fontSizeLabels = 10
fontSizeTicks = 9
fontSizePanel = 12
fontSizeTitles = 12
fontSizeNS = 10
fontSizeStars = 9

# Significance star placement
starHeightFactor = 0.2
starGapFactor = 0.3
starYfactor = 0.1

dotEdgeColor = '0.5'

rasterMarkerSize = 3 # Raster maerker size

colors = {}
#colors['blueLaser'] = cp.TangoPalette['SkyBlue1']
colors['sound'] = cp.TangoPalette['Butter2']
colors['ft'] = cp.TangoPalette['ScarletRed1']
colors['vot'] = cp.TangoPalette['SkyBlue2']
colors['audP'] = cp.TangoPalette['SkyBlue1']
colors['audD'] = cp.TangoPalette['Orange1']
colors['audV'] = cp.TangoPalette['Aluminium4']
