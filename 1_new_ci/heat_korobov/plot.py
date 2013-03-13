#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 11,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [5,4]}
py.rcParams.update(params)
import numpy as np

kappa = 2.129475
w=2

files=('R10','R16','R20','R26','R40')

for f in files:
	data = py.loadtxt(f+'.txt')
	x,y,error = data.T
	py.errorbar(x, (y-kappa)/kappa, yerr=error/kappa, lw=w, elinewidth=w-1, label=r'$\mathrm{'+f+'}$')

py.axhline(color='black')
py.xlabel(r'$N_\mathrm{kor}\times 10^5$')
py.ylabel(r'$\displaystyle\frac{\lambda-\lambda_\mathrm{ref}}{\lambda_\mathrm{ref}}$')
py.xlim(1.5,90)
py.legend(loc='lower left')
py.semilogx()
#py.loglog()

py.tight_layout()
py.savefig('convergence.pdf')
#py.show()
