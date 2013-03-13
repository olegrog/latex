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

files=('old','base','b2','Nx20-b2','minE')

for f in files:
	data = py.loadtxt(f+'.txt')
	x,y,error = data.T
        py.errorbar(x, (y-kappa)/kappa, yerr=error/kappa, lw=w, elinewidth=w-1, label=r'$\mathrm{'+f+'}$')

py.axhline(color='black')
py.xlabel(r'$N_R$')
py.ylabel(r'$\displaystyle\frac{\lambda-\lambda_\mathrm{ref}}{\lambda_\mathrm{ref}}$')
py.xlim(9,62)
py.legend((r'$\mathrm{old\ CI\ module}$', r'$\mathrm{base}$', r'$\int\dots db^2$', r'$N_x=20,db^2$', r'$\min(E)$'))
#py.semilogx()
#py.loglog()

py.tight_layout()
py.savefig('convergence.pdf')
#py.show()
