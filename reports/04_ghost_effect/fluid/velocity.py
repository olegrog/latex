#!/usr/bin/env python

import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

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

x,y = np.loadtxt('../transport/Y1_2.txt').T
Y1 = interp1d(x, y, kind='cubic')

L   = 0.5
tau = 0.5
Kn  = float(sys.argv[1])
k   = np.sqrt(np.pi)/2*Kn

def load_data():
    U,V,X,Y = data[2], data[3], data[5], data[6]
    N = int(np.sqrt(X.size))
    X = np.reshape(X,(N,N))
    Y = np.flipud(L-np.reshape(Y,(N,N)))
    x,y = X[0,:], Y[:,0]
    U = np.flipud(np.reshape(U,(N,N)))
    V = -np.flipud(np.reshape(V,(N,N)))
    if k:
        T=1-tau*np.cos(2*np.pi*x)
        dTdx=2*np.pi*tau*np.sin(2*np.pi*x)
        U -= np.sqrt(T)*dTdx*Y1(np.outer(y,1/T)/k)
    return U,V,x,y

data = np.loadtxt('asymptotic_N100.txt').T
U,V,x,y = load_data()

magU = np.sqrt(U*U+V*V)
cmap = plt.cm.get_cmap('coolwarm')
CF = plt.contourf(x, y, magU, cmap=cmap)
plt.colorbar(CF, use_gridspec=True)

# streamplot version (need matplotlib 1.2.1, use only evenly grid)
data = np.loadtxt('asymptotic_N50.txt').T
U,V,x,y = load_data()
plt.streamplot(x, y, U, V, color='k', density=.8, minlength=.4)

# quiver version
#k=5
#s=2
#Q = plt.quiver(x[s::k], x[s::k], U[s::k,s::k], V[s::k,s::k], scale=15)
#qk = plt.quiverkey(Q, -0.07, .7, 1, r'$1$', labelpos='N', labelsep=0.05,
#        fontproperties={'weight': 'bold', 'size': '8'})

plt.text(-.07, .4, r'$u_{i1}$', fontsize=11)
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
