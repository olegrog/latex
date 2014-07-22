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

from scipy.interpolate import interp1d
x,y = py.loadtxt('../transport/Y1_2.txt').T
Y1 = interp1d(x, y, kind='cubic')

L   = 	0.5
tau = 	0.5
Kn  = 	0.01
k   = 	np.sqrt(np.pi)/2*Kn

def load_data():
	U,V,X,Y = data[2], data[3], data[5], data[6]
	N = np.sqrt(X.size)
	X = np.reshape(X,(N,N))
	Y = np.flipud(L-np.reshape(Y,(N,N)))
	x,y = X[0,:], Y[:,0]
	T=1-tau*np.cos(2*np.pi*x)
	dTdx=2*np.pi*tau*np.sin(2*np.pi*x)
	U = np.flipud(np.reshape(U,(N,N)))
	U -= np.sqrt(T)*dTdx*Y1(np.outer(y,1/T)/k)
	V = -np.flipud(np.reshape(V,(N,N)))
	return U,V,x,y

data = py.loadtxt('asymptotic_N100.txt').T
U,V,x,y = load_data()

magU = np.sqrt(U*U+V*V)
cmap = py.cm.get_cmap('coolwarm')
CF = py.contourf(x, y, magU, cmap=cmap)
py.colorbar(CF, use_gridspec=True)

# streamplot version (need matplotlib 1.2.1, use only evenly grid)
data = py.loadtxt('asymptotic_N50.txt').T
U,V,x,y = load_data()
py.streamplot(x, y, U, V, color='k', density=.8, minlength=.4)

# quiver version
#k=5
#s=2
#Q = py.quiver(x[s::k], x[s::k], U[s::k,s::k], V[s::k,s::k], scale=15)
#qk = py.quiverkey(Q, -0.07, .7, 1, r'$1$', labelpos='N', labelsep=0.05,
#		fontproperties={'weight': 'bold', 'size': '8'})

py.text(-.07, .4, r'$u_{i1}$', fontsize=11)
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
