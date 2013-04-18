#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [4,3]}
py.rcParams.update(params)
import numpy as np

N=50
L=0.5
x = np.arange(0+L/2/N,L,L/N)
U, V = py.loadtxt('U_as.txt').T
U = np.flipud(np.reshape(U,(50,50)))
V = -np.flipud(np.reshape(V,(50,50)))
magU = np.sqrt(U*U+V*V)

cmap = py.cm.get_cmap('coolwarm')
CF = py.contourf(x, x, magU, cmap=cmap)
py.colorbar(CF, use_gridspec=True)

# streamplot version (need matplotlib 1.2.1)
py.streamplot(x, x, U, V, color='k', density=.8)

# quiver version
#k=5
#s=2
#Q = py.quiver(x[s::k], x[s::k], U[s::k,s::k], V[s::k,s::k], scale=15)
#qk = py.quiverkey(Q, -0.07, .7, 1, r'$1$', labelpos='N', labelsep=0.05,
#		fontproperties={'weight': 'bold', 'size': '8'})

py.text(-.05, .4, r'$u_i$', fontsize=11)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$z$', labelpad=-5, rotation=0)
py.xlim(0,0.5)
py.ylim(0,0.5)
ax = py.axes()
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
py.tick_params(axis='both', direction='out')
py.savefig('velocity.pdf', bbox_inches='tight')
