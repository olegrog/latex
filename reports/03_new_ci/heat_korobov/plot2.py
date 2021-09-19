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

files=('R10','R14','R16','R20','R26','R40')

for f in files:
	data = np.loadtxt(f+'.txt')
	x,y,error = data.T
	plt.plot(x, error*np.sqrt(41), lw=w, label=r'$\mathrm{'+f+'}$')

plt.axhline(color='black')
plt.xlabel(r'$N_\mathrm{kor}\times 10^5$')
plt.ylabel(r'$\sigma(\lambda)$')
plt.xlim(1.5,90)
plt.legend(loc='upper right')
plt.semilogx()
#plt.loglog()

plt.tight_layout()
plt.savefig('korobov.pdf')
#plt.show()
