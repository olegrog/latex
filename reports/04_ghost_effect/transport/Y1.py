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
    'figure.figsize': [3,3]
}
plt.rcParams.update(params)

data = np.loadtxt('Y1_2.txt')
x,y = data.T
p, = plt.plot(x, y, '-', lw=2)

plt.xlabel(r'$\eta$', labelpad=-1)
plt.semilogy()
plt.axes().set_xticks([0,5,10,15])
plt.axes().set_yticks([1,1e-2,1e-4])
plt.ylabel(r'$\displaystyle\frac{Y_1}2$', y=.8, labelpad=-3, rotation=0)
plt.xlim(0,15)

plt.tight_layout()
plt.savefig('Y1.pdf')
