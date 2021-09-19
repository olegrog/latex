#!/usr/bin/env python

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.use('pgf')
params = {
    'axes.labelsize': 11,
    'font.size': 11,
    'legend.fontsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'pgf.rcfonts': False,
    'figure.figsize': [5,4]
}
plt.rcParams.update(params)

kappa = 2.129475
w=2

files=('old','base','b2','Nx20-b2','minE')

for f in files:
    data = np.loadtxt(f+'.txt')
    x,y,error = data.T
    plt.errorbar(x, (y-kappa)/kappa, yerr=error/kappa, lw=w, elinewidth=w-1, label=r'$\mathrm{'+f+'}$')

plt.axhline(color='black')
plt.xlabel(r'$N_R$')
plt.ylabel(r'$\displaystyle\frac{\lambda-\lambda_\mathrm{ref}}{\lambda_\mathrm{ref}}$')
plt.xlim(9,62)
plt.legend((r'$\mathrm{old\ CI\ module}$', r'$\mathrm{base}$', r'$\int\dots db^2$', r'$N_x=20,db^2$', r'$\min(E)$'))
#plt.semilogx()
#plt.loglog()

plt.tight_layout()
plt.savefig('convergence.pdf')
#plt.show()
