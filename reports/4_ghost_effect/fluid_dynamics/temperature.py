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

L = 0.5

#data = py.loadtxt('heat.txt').T
#T,X,Y = data[0], data[1], data[2]
data = py.loadtxt('asymptotic_N50.txt').T
T,X,Y = data[1], data[5], data[6]
N = np.sqrt(X.size)
X = np.reshape(X,(N,N))
Y = np.flipud(L-np.reshape(Y,(N,N)))
x,y = X[0,:], Y[:,0]
T = np.flipud(np.reshape(T,(N,N)))

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
py.savefig('temperature.pdf', bbox_inches='tight')
