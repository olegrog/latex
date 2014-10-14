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
import scipy.optimize
import functools
import sys

data = py.loadtxt(sys.argv[1] + '.txt')
alpha, maxU = data.T[0:2]
tau = alpha - 1
py.plot(tau, maxU, '-', lw=1.5)

def func(n, x, A):
    return A*x**n;

def plot_dash_line(n, ind):
    X = np.array([ tau.min(), tau.max() ])
    Y = func(n, X, np.exp( np.log(maxU[ind]) - n*np.log(tau[ind]) ))
    py.plot(X, Y, 'k--', lw=.5, label=r'$\phi = \pi/2$')

plot_dash_line(3, 0)
plot_dash_line(1.5, -1)

py.xlabel(r'$\tau$', labelpad=5)
py.ylabel(r'$(|u_{i1}|)_{\mathrm{max}}$', y=.8, labelpad=0, rotation=0)
py.loglog()
py.xlim(1e-2, 1e2)
py.ylim(1e-8, 1e4)
ax = py.axes()

from matplotlib.ticker import LogLocator
ax.yaxis.set_major_locator(LogLocator(1e4))

py.savefig('U-tau-' + sys.argv[1] + '.pdf', bbox_inches='tight')

