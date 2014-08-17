#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np

def plot(filename, factor, fmt, lw, label):
    d,Q,F = py.loadtxt(filename).T
    F = F * factor / 1e3
    py.plot(d, F, fmt, lw=lw, label=label)

plot('cylinders.txt', -2, 'b-', 1, r'$\mathrm{cylinders}$')
plot('spheres.txt', -100, 'g-', 1, r'$\mathrm{spheres}$')

py.xlabel(r'$d$', labelpad=-3)
py.ylabel(r'$-\displaystyle\frac{F_{i2}d_i}{d}\times10^{-3}$', labelpad=0)
py.xlim(0, 0.9)
py.ylim(0, 6)
ax = py.axes()
ax.set_xticks([0,0.3,0.6,0.9])
ax.set_yticks([0,2,4,6])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.yaxis.set_minor_locator(MultipleLocator(1))

py.savefig('force.pdf', bbox_inches='tight')

