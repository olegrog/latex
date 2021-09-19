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
    'figure.figsize': [3.6,2.8]
}
plt.rcParams.update(params)

L  = 0.5
Kn = float(sys.argv[1])
k  = np.sqrt(np.pi)/2*Kn

def load_data():
    rho,U,V,X,Y = data[1], data[2], data[3], data[4], data[5]
    N = int(np.sqrt(X.size))
    M = X.max()
    X = np.reshape(X,(N,N)) * L/M
    Y = np.reshape(Y,(N,N)) * L/M
    x,y = X[0,:], Y[:,0]
    U = np.reshape(U/rho,(N,N)) / k
    V = np.reshape(V/rho,(N,N)) / k
    return U,V,x,y

data = np.loadtxt(sys.argv[2]).T
U,V,x,y = load_data()

magU = np.sqrt(U*U+V*V)
cmap = plt.cm.get_cmap('coolwarm')
CF = plt.contourf(x, y, magU, cmap=cmap)
plt.colorbar(CF, use_gridspec=True)

# streamplot version (need matplotlib 1.2.1, use only evenly grid)
plt.streamplot(x, y, U, V, color='k', density=.8, minlength=.4)

plt.text(-.06, .37, r'$\displaystyle\frac{u_i}k$', fontsize=11)
plt.xlabel(r'$x$', labelpad=-5)
plt.ylabel(r'$z$', labelpad=-5, rotation=0)
plt.xlim(0,0.5)
plt.ylim(0,0.5)
ax = plt.axes()
ax.set_aspect('equal')
ax.set_yticks([0,0.5])
ax.set_xticks([0,0.5])
plt.tick_params(axis='both', direction='out')
plt.savefig(sys.argv[3], bbox_inches='tight')
