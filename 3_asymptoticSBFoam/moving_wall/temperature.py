#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [6,5]}
py.rcParams.update(params)
import numpy as np

L = 0.5

#data = py.loadtxt('heat.txt').T
#T,X,Y = data[0], data[1], data[2]
data = py.loadtxt('asym.txt').T
T,X,Y = data[1], data[5], data[6]
N = 50
dim = (N+1,2*N+1)
X = np.reshape(X,dim)
Y = np.reshape(Y,dim)
T = np.reshape(T,dim)
x,y = X[0,:], Y[:,0]

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

py.text(-.05, .4, r'$T$', fontsize=10)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$z$', labelpad=-5, rotation=0)
py.xlim(0,1)
py.ylim(0,0.5)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,1])
py.tick_params(axis='both', direction='out')
py.savefig('T.pdf', bbox_inches='tight')
