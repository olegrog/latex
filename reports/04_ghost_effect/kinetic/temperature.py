#!/usr/bin/env python

import sys
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

L = 0.5

data = np.loadtxt(sys.argv[1]).T
T,X,Y = data[0], data[4], data[5]
N = int(np.sqrt(X.size))
M = X.max()
X = np.reshape(X,(N,N)) * L/M
Y = np.reshape(Y,(N,N)) * L/M
x,y = X[0,:], Y[:,0]
T = np.reshape(T,(N,N))

levels=np.arange(.5,1.6,.1)
cmap = plt.cm.get_cmap('coolwarm')
plt.contourf(x, y, T, levels=levels, cmap=cmap)

levels = np.append(levels, 1.05)
levels.sort()
CS = plt.contour(x, y, T, levels=levels, colors='k', linewidths=1)
plt.clabel(CS, levels[1::1],
          inline=True,
          use_clabeltext=True,
          fmt='%g',
          fontsize=8)

plt.text(-.05, .4, r'$T$', fontsize=11)
plt.xlabel(r'$x$', labelpad=-5)
plt.ylabel(r'$z$', labelpad=-5, rotation=0)
plt.xlim(0,0.5)
plt.ylim(0,0.5)
ax = plt.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
plt.tick_params(axis='both', direction='out')
plt.savefig(sys.argv[2], bbox_inches='tight')
