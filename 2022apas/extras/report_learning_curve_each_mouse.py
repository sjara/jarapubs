"""
Plot the learning curve for each mouse.
"""

import os
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from jaratoolbox import settings
from jaratoolbox import extraplots
from jaratoolbox import extrastats
import studyparams
import studyutils
import figparams
from importlib import reload
reload(figparams)
reload(studyutils)



dframe = studyutils.load_stage3(excludeAntibias=False)
