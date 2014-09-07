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

def pow1(x,A):
    return A*x**1

def plot(filename, factor, fmt, lw, label):
    d, Q, F = py.loadtxt(filename).T
    F = F * factor
    py.plot(d, F, fmt, lw=lw, label=label)
    ran = np.where(d < .05)
    p, pcov = scipy.optimize.curve_fit(pow1, d[ran], F[ran], p0=[1])
    print ran, p
    X = np.logspace(np.log10(1e-2), np.log10(0.6), num=50)
    Y = pow1(X, p[0])
    py.plot(X, Y, 'k--', lw=.5,)

plot('cylinders.txt', -2, '-', 1.5, r'$\mathrm{cylinders}\;\alpha=5$')
plot('spheres.txt', -100, '-', 1.5, r'$\mathrm{spheres}\;\alpha=5$')

py.xlabel(r'$d$', labelpad=-5)
py.ylabel(r'$-\displaystyle\oint p_0 F_{x2} \mathrm{d}S$', labelpad=0)
py.semilogx()
py.semilogy()
ax = py.axes()
ax.set_xticks([1e-2, 1e-1, 1])
ax.set_xticklabels([r'$10^{-2}$', '', r'$10^0$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

py.savefig('forces.pdf', bbox_inches='tight')

