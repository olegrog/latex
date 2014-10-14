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
import scipy.optimize
import functools

def pow(n, x, A):
    return A*x**n

def plot(filename, factor, label, n):
    d, Q, F, eps = py.loadtxt(filename).T
    F = F * factor
    d, F = 1-d, 1/F
    py.plot(d, F, '-', lw=1.5, label=label)
    func = functools.partial(pow, n)
    X = np.logspace(np.log10(1e-2), np.log10(1), num=50)
    Y = func(X, F.min()/func(d.min(), 1))
    py.plot(X, Y, 'k--', lw=.5,)

plot('cylinders.txt', 2, r'$\mathrm{cylinders}\;\tau=4$', 1.5)
plot('spheres.txt', 100, r'$\mathrm{spheres}\;\tau=4$', 1)

py.xlabel(r'$1-d$', labelpad=-5)
py.ylabel(r'$\displaystyle \left(p_0\oint_S F_{x2} \mathrm{d}S \right)^{-1}$', labelpad=0)
py.loglog()
py.ylim(1e-5, 1e-1)
ax = py.axes()
ax.set_xticks([1e-2, 1e-1, 1])
ax.set_xticklabels([r'$10^{-2}$', '', r'$10^0$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

py.savefig('forces-close.pdf', bbox_inches='tight')

