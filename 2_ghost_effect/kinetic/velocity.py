#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 11,
          'text.fontsize': 11,
          'legend.fontsize': 12,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'text.usetex': True,
          'figure.figsize': [3.6,2.8]}
py.rcParams.update(params)
import numpy as np

L  = 0.5
Kn = 0.01
k  = np.sqrt(np.pi)/2*Kn

def load_data():
	rho,U,V,X,Y = data[1], data[2], data[3], data[4], data[5]
	N = np.sqrt(X.size)
	M = X.max()
	X = np.reshape(X,(N,N)) * L/M
	Y = np.reshape(Y,(N,N)) * L/M
	x,y = X[0,:], Y[:,0]
	U = np.reshape(U/rho,(N,N)) / k 
	V = np.reshape(V/rho,(N,N)) / k
	return U,V,x,y

data = py.loadtxt('kn001.txt').T
U,V,x,y = load_data()

magU = np.sqrt(U*U+V*V)
cmap = py.cm.get_cmap('coolwarm')
CF = py.contourf(x, y, magU, cmap=cmap)
py.colorbar(CF, use_gridspec=True)

# streamplot version (need matplotlib 1.2.1, use only evenly grid)
py.streamplot(x, y, U, V, color='k', density=.8, minlength=.4)

py.text(-.06, .37, r'$\displaystyle\frac{u_i}k$', fontsize=11)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$z$', labelpad=-5, rotation=0)
py.xlim(0,0.5)
py.ylim(0,0.5)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
py.tick_params(axis='both', direction='out')
py.savefig('velocity.pdf', bbox_inches='tight')
