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
import sys

def func(n, x, A):
    return A*x**n;

def plot_dash_line(n, ind):
    X = np.array([ x.min(), x.max() ])
    Y = func(n, X, np.exp( np.log(F[ind]) - n*np.log(x[ind]) ))
    py.plot(X, Y, 'k--', lw=.5)

def plot(name, label, column, factor):
    data = py.loadtxt(name + '.txt')
    alpha, F = data.T[0:column:column-1]
    x, F = alpha-1, factor*F
    py.plot(x, F, '-', lw=1, label=label)
    return x, F

if sys.argv[2] == 'inner':
    factor = 2
    column = 4
    signP = -1
if sys.argv[2] == 'outer':
    factor = -2
    column = 8
    signP = 1

x, F = plot(sys.argv[1], r'$\mathrm{term\;}' + ("% d" % signP)[0] + r'p_2^\dag$', column+1, signP*factor)
if sys.argv[2] == 'outer':
    plot_dash_line(3, 0)
x, F = plot(sys.argv[1], r'$\mathrm{term\;}-\partial u_{i1}/\partial x_j$', column+2, -factor)
plot_dash_line(3, 0)
plot_dash_line(2, -1)
x, F = plot(sys.argv[1], r'$\mathrm{term\;}(\partial T_0/\partial x_k)^2$', column+3, factor)
x, F = plot(sys.argv[1], r'$\mathrm{total\;value}$', column, factor)

py.xlabel(r'$\tau$', labelpad=5)
py.ylabel(r'$' + ("% d" % factor)[0] + r'p_0\displaystyle\oint_S F_{x2}\mathrm{d}S$', labelpad=-5)
py.loglog()
py.ylim(1e-6, 1e6)
ax = py.axes()
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

from matplotlib.ticker import LogLocator
ax.yaxis.set_major_locator(LogLocator(1e2))

py.savefig('F-tau-' + sys.argv[1] + '-' + sys.argv[2] + '.pdf', bbox_inches='tight', transparent=True)

