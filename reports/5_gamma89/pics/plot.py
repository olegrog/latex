#!/usr/bin/env python
# Usage: ./plot.py <file> <column> <ylabel> <output>

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,2]}
py.rcParams.update(params)
import numpy as np
import sys
import math

corrector = lambda x: x
kind = 'simple'
y_coord = .8

def plot(filename, column, label):
    data = py.loadtxt(filename).T
    X, Y = data[0], data[column]
    mask = (X >= xmin) * (X <= xmax)
    X, Y = X[mask], corrector(Y[mask])
    py.ylabel(r'$' + label + r'$', y=y_coord, labelpad=8, rotation=0)
    py.plot(X, Y, '-', lw=1)
    return np.min(Y), np.max(Y)

correctors = {
    'log': lambda y: np.abs(y)
}

if len(sys.argv) > 5:
    kind = sys.argv[5]
    corrector = correctors[kind]

if len(sys.argv) > 6:
    y_coord = float(sys.argv[6])

xmin, xmax = 0.0, 5.0
ymin, ymax = plot(sys.argv[1], sys.argv[2], sys.argv[3])

py.xlabel(r'$\zeta$', labelpad=-5)
py.xlim(xmin, xmax)
ax = py.axes()
ax.axhline(lw=.5, c='k', ls=':')

L = max(0, ymax) - min(0, ymin)
tick = 10**math.floor(np.log10(L))
if L/tick > 6:
    tick *= 2
elif L/tick < 1.5:
    tick /= 5
elif L/tick < 3:
    tick /= 2

if kind == 'log':
    py.semilogy()
else:
    from matplotlib.ticker import MultipleLocator
    ax.yaxis.set_major_locator(MultipleLocator(tick))

py.savefig(sys.argv[4], bbox_inches='tight')

