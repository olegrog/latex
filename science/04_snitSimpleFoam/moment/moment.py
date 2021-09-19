#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np

def plot(filename, factor, fmt, lw, label):
    beta, Q, M = py.loadtxt(filename).T
    beta = .5 - beta
    M = M * factor
    py.plot(beta, M, fmt, lw=lw, label=label)

plot('elliptic.txt', 2, '-', 1.5, r'$\alpha = 5$')
#plot('inverse.txt', 2, '-', 1.5, r'$\alpha = -5$')

py.xlabel(r'$\beta$', labelpad=-5)
py.ylabel(r'$\displaystyle p_0 \oint M_{z2} \mathrm{d}S$', labelpad=0)
ax = py.axes()
ax.set_xticks([0, 0.5])
ax.set_xticklabels([r'$0$', r'$\pi/2$'])
handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=0)
ax.axhline(lw=.5, c='k', ls=':')

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(0.25))
#ax.yaxis.set_minor_locator(MultipleLocator(10))

py.savefig('moment-beta.pdf', bbox_inches='tight', transparent=True)

