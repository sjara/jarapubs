import numpy as np
import matplotlib.pyplot as plt

def scatter_hist(x, y, ax, ax_histx, ax_histy, binwidth = None):
    ...:     if binwidth is None:
    ...:         binwidth = 0.1
    ...:     xymax = max(np.max(np.abs(x)), np.max(np.abs(y)))
    ...:     lim = (int(xymax/binwidth) + 1) * binwidth
    ...:     bins = np.arange(-lim, lim+binwidth, binwidth)
    ...:     ax_histx.tick_params(axis = "x", labelbottom = True)
    ...:     ax_histy.tick_params(axis = "y", labelleft = True)
    ...:     ax.scatter(x, y)
    ...:     ax_histx.hist(x, bins = bins)
    ...:     ax_histy.hist(y, bins = bins, orientation = 'horizontal')
