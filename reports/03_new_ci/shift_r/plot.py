#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [5,4]}
py.rcParams.update(params)
import numpy as np

kappa = 2.129475
w=2

data = py.loadtxt('R16_k5.txt')
x,y1,y2 = data.T
py.plot(x, (y1-kappa)/kappa, lw=w, label=r'$\mathrm{normal}$')
py.plot(x, (y2-kappa)/kappa, lw=w, label=r'$\mathrm{shift}$')

py.axhline(color='black')
py.xlabel(r'$time$')
py.ylabel(r'$\displaystyle\frac{\lambda-\lambda_\mathrm{ref}}{\lambda_\mathrm{ref}}$')
#py.xlim(9,62)
py.legend(loc='upper left')
#py.semilogx()
#py.loglog()

py.tight_layout()
py.savefig('convergence.pdf')
#py.show()
