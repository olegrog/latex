#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np
import sys

L = 0.5

data = py.loadtxt(sys.argv[1]).T
T,X,Y = data[0], data[4], data[5]
N = np.sqrt(X.size)
M = X.max()
X = np.reshape(X,(N,N)) * L/M
Y = np.reshape(Y,(N,N)) * L/M
x,y = X[0,:], Y[:,0]
T = np.reshape(T,(N,N))

levels=np.arange(.5,1.6,.1)
cmap = py.cm.get_cmap('coolwarm')
py.contourf(x, y, T, levels=levels, cmap=cmap)

levels = np.append(levels, 1.05)
CS = py.contour(x, y, T, levels=levels, colors='k', linewidths=1)
py.clabel(CS, levels[1::1],
          inline=True,
          use_clabeltext=True,
          fmt='%g',
          fontsize=8)

py.text(-.05, .4, r'$T$', fontsize=11)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$z$', labelpad=-5, rotation=0)
py.xlim(0,0.5)
py.ylim(0,0.5)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
py.tick_params(axis='both', direction='out')
py.savefig(sys.argv[2], bbox_inches='tight')
