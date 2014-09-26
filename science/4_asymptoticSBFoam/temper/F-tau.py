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

def func(n, x, A):
    return A*x**n;

def plot_dash_line(n, ind):
    X = np.array([ x.min(), x.max() ])
    Y = func(n, X, np.exp( np.log(F[ind]) - n*np.log(x[ind]) ))
    py.plot(X, Y, 'k--', lw=.5)

def plot(filename, label, factor):
    data = py.loadtxt(filename)
    alpha, F = data.T[0:4:3]
    x, F = alpha-1, factor*F
    py.plot(x, F, '-', lw=1.5, label=label)
    return x, F

x, F = plot('cylinders.txt', r'$\mathrm{cylinders}$', -2)
plot_dash_line(2, 0)
plot_dash_line(3, -1)

#x, F = plot('spheres.txt', r'$\mathrm{spheres}$', -100)
#plot_dash_line(1, 0)
#plot_dash_line(3, -1)

py.xlabel(r'$\tau$', labelpad=5)
py.ylabel(r'$-p_0F_{x2}$', y=.8, labelpad=0, rotation=0)
py.loglog()
py.xlim(1e-2, 1e2)
py.ylim(1e-4, 1e8)
ax = py.axes()
#handles, labels = ax.get_legend_handles_labels()
#ax.legend(handles, labels, loc=0)

from matplotlib.ticker import LogLocator
ax.yaxis.set_major_locator(LogLocator(1e4))

py.savefig('F-tau-noncoaxial.pdf', bbox_inches='tight')

