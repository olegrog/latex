#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 9,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [5,3]}
py.rcParams.update(params)
import numpy as np

def plot(filename, label):
    x,y = py.loadtxt(filename).T
    py.plot(x, y, '-', lw=1, label=label)

plot('p_0',  r'$p^\dag_2$')
plot('Ux_0', r'$u_{x1}$')
plot('Uy_0', r'$u_{y1}$')
plot('T_0',  r'$T_0$')

py.xlabel(r'$n$', labelpad=5)
py.ylabel(r'$\varepsilon$', y=.75, labelpad=5, rotation=0)
ax = py.axes()
ax.semilogy()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

from matplotlib.ticker import LogLocator
ax.yaxis.set_major_locator(LogLocator(1e2))

py.savefig('residual.pdf', bbox_inches='tight', transparent=True)

