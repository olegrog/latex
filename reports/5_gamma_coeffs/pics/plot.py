#!/usr/bin/env python
# Usage: ./plot.py <file> <column> <ylabel> <output> [<corrector>] [<y_coord>] [<box>]

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
from matplotlib.ticker import MultipleLocator, MaxNLocator, LinearLocator

xmin, xmax = 0.0, 5.0
corrector = lambda y: y
kind = 'simple'
y_coord = .8
inside = False

def plot(filename, column, label):
    data = py.loadtxt(filename).T
    X, Y = data[0], data[column]
    mask = (X >= xmin) * (X <= xmax)
    X, Y = X[mask], corrector(Y[mask])
    py.ylabel(r'$' + label + r'$', y=y_coord, labelpad=8, rotation=0)
    py.plot(X, Y, '-', lw=1)
    py.xlabel(r'$\zeta$', labelpad=-5)
    py.xlim(xmin, xmax)
    ax = py.axes()
    ax.axhline(lw=.5, c='k', ls=':')
    specify_tics(np.min(Y), np.max(Y))
    if inside:
        print "Plot inside plot"
        ax = py.axes([.25, .45, .4, .4])
        mask = X <= float(sys.argv[7])
        py.plot(X[mask], Y[mask], '-', lw=1)
        ax.tick_params(axis='both', which='major', labelsize=8)
        ax.axhline(lw=.5, c='k', ls=':')
        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(LinearLocator(3))
        ymin, _, ymax = ax.get_yticks()
        if (ymax + ymin) < (ymax-ymin)/5:
            y = max(ymax, -ymin)
            py.ylim(-y, y)

def specify_tics(ymin, ymax):
    L = max(0, ymax) - min(0, ymin)
    tick = 10**math.floor(np.log10(L))
    print L, tick, L/tick
    if L/tick > 6:
        tick *= 2
    elif L/tick < 2:
        tick /= 2
    elif L/tick < 1.2:
        tick /= 5
    elif L/tick < 3:
        tick /= 2
    print tick

    if kind == 'log':
        py.semilogy()
    else:
        ax = py.axes()
        ax.yaxis.set_major_locator(MultipleLocator(tick))


correctors = {
    'log': lambda y: np.abs(y),
    'none': lambda y: y
}

if len(sys.argv) > 5:
    kind = sys.argv[5]
    corrector = correctors[kind]

if len(sys.argv) > 6:
    if sys.argv[6] != 'default':
        y_coord = float(sys.argv[6])

if len(sys.argv) > 7:
    inside = True

plot(sys.argv[1], sys.argv[2], sys.argv[3])

py.savefig(sys.argv[4], bbox_inches='tight')

