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

def plot(filename, factor, fmt, lw, label):
    d,Q,F = py.loadtxt(filename).T
    F = F * factor
    py.plot(d, F, fmt, lw=lw, label=label)

plot('cylinders.txt', -2, '-', 1, r'$\mathrm{cylinders}\;\alpha=5$')
plot('spheres.txt', -100, '-', 1, r'$\mathrm{spheres}\;\alpha=5$')

py.xlabel(r'$d$', labelpad=-5)
py.ylabel(r'$-\displaystyle\oint F_{x2} \mathrm{d}S$', labelpad=0)
py.semilogx()
py.semilogy()
ax = py.axes()
ax.set_xticks([1e-2, 1e-1, 1])
ax.set_xticklabels([r'$10^{-2}$', '', r'$10^0$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

py.savefig('forces.pdf', bbox_inches='tight')

