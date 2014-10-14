#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np
import mpmath

major, minor = .7, .3
L = 2*major*mpmath.ellipe(1-np.square(minor/major))
xmin, xmax = 0, 1
def plot(filename, fmt, lw, label, factor=1):
    x,y = py.loadtxt(filename).T
    py.plot(x/L, y, fmt, lw=lw, label=label)

plot('beta500.raw', '-', 1, r'$\beta = 0$', -1)
plot('beta375.raw', '-', 1, r'$\beta = \pi/8$', -1)
plot('beta250.raw', '-', 1, r'$\beta = \pi/4$', -1)
plot('beta125.raw', '-', 1, r'$\beta = 3\pi/8$', -1)
plot('beta000.raw', '-', 1, r'$\beta = \pi/2$', -1)

py.xlabel(r'$s$', labelpad=-5)
py.ylabel(r'$p_0 M_{z2}$', y=.75, labelpad=0, rotation=0)
py.xlim(xmin, xmax)
ax = py.axes()
ax.set_xticks([xmin, xmax])
#ax.set_yticks([-400, 0, 400])
ax.set_xticklabels([r'$0$', r'$L_\mathrm{ell}/2$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)
ax.axhline(lw=.5, c='k', ls=':')

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(0.5*xmax))
ax.yaxis.set_minor_locator(MultipleLocator(100))

py.savefig('profiles.pdf', bbox_inches='tight')

