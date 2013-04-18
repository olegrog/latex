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

N=50
L=0.5
x = np.arange(0+L/2/N,L,L/N)
T = np.flipud(np.reshape(py.loadtxt('T_he.txt'),(50,50)))

levels=np.arange(.5,1.6,.1)
cmap = py.cm.get_cmap('coolwarm')
py.contourf(x, x, T, levels=levels, cmap=cmap)

levels = np.append(levels, 1.05)
CS = py.contour(x, x, T, levels=levels, colors='k', linewidths=1)
py.clabel(CS, levels[1::1],
          inline=True,
	  use_clabeltext=True,
          fmt='%g',
          fontsize=8)

py.text(-.05, .4, r'$T$', fontsize=10)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$z$', labelpad=-5, rotation=0)
py.xlim(0,0.5)
py.ylim(0,0.5)
ax = py.axes()
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
py.tick_params(axis='both', direction='out')
py.tight_layout()
py.savefig('temperature.pdf')
