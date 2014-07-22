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

data = py.loadtxt('alpha.txt')
x,y = data.T
xmin = x.min()
xmax = x.max()

from scipy.optimize import curve_fit
def func(x,A,k):
	return A*(x-1)**3

ran = np.where((x-1)<.05)
p, pcov = curve_fit(func, x[ran], y[ran], p0=[0.02,3])
print ran
print p

X = np.logspace(np.log10(xmin), np.log10(10), num=50)
Y = func(X, p[0], p[1])

py.plot(x-1, y, '-', lw=2)
py.plot(X-1, Y, '--', lw=1)

py.xlabel(r'$\alpha$', labelpad=-3)
py.ylabel(r'$(u_{i1})_{\mathrm{max}}$', y=.8, labelpad=-15, rotation=0)
py.xlim(xmin-1, xmax-1)
ax = py.axes()
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks([1e-1,1e+1])
ax.set_yticks([1e-9,1e-3,1e+2])

py.savefig('alpha.pdf', bbox_inches='tight')

