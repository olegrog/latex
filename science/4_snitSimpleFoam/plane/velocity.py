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
import sys

from scipy.interpolate import interp1d
x,y = py.loadtxt('../transport/Y1_2.txt').T
Y1 = interp1d(x, y, kind='cubic')

L     = 0.5
alpha = 0.5
beta  = 2./np.sqrt(np.pi)
Kn    = float(sys.argv[1])
k     = np.sqrt(np.pi)/2*Kn

data = py.loadtxt('asym.raw').T
X0, Y0, U0, V0 = data[0], data[1], data[3], data[4]
N = 50
dim = (N+1, 2*N+1)
X = np.reshape(X0, dim)
Y = np.reshape(Y0, dim)
x,y = X[0,:], Y[:,0]
T = 1-alpha*np.cos(2*np.pi*x)
dTdx = 2*np.pi*alpha*np.sin(2*np.pi*x)
U = np.reshape(U0,dim)
if k:
    U -= np.sqrt(T)*dTdx*Y1(np.outer(y,1/T)/k)
U0 = np.reshape(U, ((2*N+1)*(N+1),))
V = np.reshape(V0,dim)

np.set_printoptions(threshold='nan')
magU = np.sqrt(U*U+V*V)
cmap = py.cm.get_cmap('coolwarm')
lmax = 3.2 # (int(magU.max()*10) + 1.)/10
lev = np.linspace(0, lmax, 11)
CF = py.contourf(x, y, magU, cmap=cmap, levels=lev)

# streamplot version (need matplotlib 1.2.1, use only evenly grid)
from matplotlib.mlab import griddata
interp = 'nn'
N = 50
xi = np.linspace(0, 1, 2*N)
yi = np.linspace(0, .5, 4*N)
U = griddata(X0, Y0, U0, xi, yi, interp=interp)
V = griddata(X0, Y0, V0, xi, yi, interp=interp)
py.streamplot(xi, yi, U, V, color='k', density=0.89, minlength=.2, arrowstyle='->')

py.text(-.07, .4, r'$u_{i1}$', fontsize=10)
py.xlabel(r'$x$', labelpad=-3)
py.ylabel(r'$z$', labelpad=-3, rotation=0)
py.xlim(0,1)
py.ylim(0,0.5)
ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,1])

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(py.axes())
cax = divider.append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax)

py.tick_params(axis='both', direction='out')
py.savefig('U' + '{0:g}'.format(Kn)[2:] + '.pdf', bbox_inches='tight')
