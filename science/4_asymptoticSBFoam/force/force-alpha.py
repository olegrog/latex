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

def pow3(x,A,k):
    return A*(x-1)**3

alpha, Q, F = py.loadtxt('alpha.txt').T
F = -2 * F
py.plot(alpha-1, F, '-', lw=1.5, label=r'$\mathrm{cylinders}\;d=0.5$')

ran = np.where((alpha-1)>10)
p, pcov = scipy.optimize.curve_fit(pow3, alpha[ran], F[ran], p0=[0.02, 3])
print ran, p
X = np.logspace(np.log10(0.1), np.log10(50), num=50)
Y = pow3(X+1, p[0], p[1])
py.plot(X, Y, 'k--', lw=.5)

py.xlabel(r'$\alpha-1$', labelpad=-5)
py.ylabel(r'$-\displaystyle\oint p_0 F_{x2} \mathrm{d}S$', labelpad=-8)
py.semilogx()
py.semilogy()
ax = py.axes()
ax.set_xticks([1e-2, 1e-1, 1, 1e1])
ax.set_yticks([1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6])
ax.set_xticklabels([r'$10^{-2}$', '', '', r'$10^1$'])
ax.set_yticklabels([r'$10^{-3}$', '', '', r'$10^0$', '', '', r'$10^3$', '', '', r'$10^6$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=0)

py.savefig('forces-alpha.pdf', bbox_inches='tight')

