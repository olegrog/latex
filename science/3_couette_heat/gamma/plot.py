#!/usr/bin/env python
# Usage: ./plot.py <file> <column> <ylabel> <output> [<corrector>] [<y_coord>] [<box>]

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 8,
          'font.size': 8,
          'legend.fontsize': 8,
          'xtick.labelsize': 8,
          'ytick.labelsize': 8,
          'text.usetex': True,
          'figure.figsize': [2,1.5]}
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
    aY, bY = np.min(Y), np.max(Y)
    pad = 7 if aY < 0 and -bY/aY < 10 else 0
    py.ylabel(r'$' + label + r'$', y=y_coord, labelpad=8-pad, rotation=0)
    py.plot(X, Y, '-', lw=1)
    py.xlabel(r'$\zeta$', labelpad=-5)
    py.xlim(xmin, xmax)
    ax = py.axes()
    ax.axhline(lw=.5, c='k', ls=':')
    specify_tics(np.min(Y), np.max(Y))
    if inside:
        print "Plot inside plot"
        ax = py.axes([.27, .47, .38, .38])
        mask = X <= float(sys.argv[7])
        py.plot(X[mask], Y[mask], '-', lw=1)
        ax.tick_params(axis='both', which='major', labelsize=7)
        ax.axhline(lw=.5, c='k', ls=':')
        ax.xaxis.set_major_locator(LinearLocator(2))
        ax.xaxis.set_minor_locator(MultipleLocator(0.5))
        ax.yaxis.set_major_locator(LinearLocator(2))
        ymin, ymax = ax.get_yticks()
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

py.savefig(sys.argv[4], bbox_inches='tight', pad_inches=0.05)

