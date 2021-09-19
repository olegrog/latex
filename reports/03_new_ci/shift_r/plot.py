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

data = np.loadtxt('R16_k5.txt')
x,y1,y2 = data.T
plt.plot(x, (y1-kappa)/kappa, lw=w, label=r'$\mathrm{normal}$')
plt.plot(x, (y2-kappa)/kappa, lw=w, label=r'$\mathrm{shift}$')

plt.axhline(color='black')
plt.xlabel(r'$\mathrm{time}$')
plt.ylabel(r'$\displaystyle\frac{\lambda-\lambda_\mathrm{ref}}{\lambda_\mathrm{ref}}$')
#plt.xlim(9,62)
plt.legend(loc='upper left')
#plt.semilogx()
#plt.loglog()

plt.tight_layout()
plt.savefig('convergence.pdf')
#plt.show()
