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

d1 = py.loadtxt('heat5.txt')
y1,x1 = d1.T
d2 = py.loadtxt('asym5.txt')
y2,x2 = d2.T
py.plot(x1, y2-y1, '-', lw=1.5)

py.xlabel(r'$x$', labelpad=-2)
py.ylabel(r'$\delta T$', y=.8, labelpad=-10, rotation=0)
py.xlim(0.3, 1.5)

ax = py.axes()
ax.set_xticks([.3, 1.5])
ax.set_yticks([-6e-3, -3e-3, 0])

from matplotlib.ticker import MultipleLocator
ml = MultipleLocator(0.2)
ax.xaxis.set_minor_locator(ml)
ml = MultipleLocator(1e-3)
ax.yaxis.set_minor_locator(ml)


py.savefig('temper.pdf', bbox_inches='tight')
