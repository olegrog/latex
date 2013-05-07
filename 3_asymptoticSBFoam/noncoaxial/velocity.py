#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 9,
          'text.fontsize': 9,
          'legend.fontsize': 9,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [6,5]}
py.rcParams.update(params)
import numpy as np

d = .5
r = 2
times = 1e3

data = py.loadtxt('data.txt').T
U,V,X,Y = data[2], data[3], data[5], data[6]
magU = np.sqrt(U*U+V*V) * times

import matplotlib.tri as tri
triang = tri.Triangulation(X,Y)

def sqr(x):
	return x*x

xmid = X[triang.triangles].mean(axis=1)
ymid = Y[triang.triangles].mean(axis=1)
mask = np.where(sqr(xmid-d) + sqr(ymid) < sqr(1), 1, 0)
triang.set_mask(mask)

#lev = np.linspace(magU.min(), magU.max(), 8)
lev = np.linspace(0, 2.6, 7)
cmap = py.cm.get_cmap('coolwarm')
CF = py.tricontourf(triang, magU, cmap=cmap, levels=lev)
#py.tricontour(X, Y, magU, levels=lev, colors='k', linewidths=.5)

from matplotlib.mlab import griddata
N = 50
interp = 'linear'
xi = np.linspace(-r, r, 2*N)
yi = np.linspace(0 , r, N  )
U = griddata(X,Y,U,xi,yi,interp=interp)
V = griddata(X,Y,V,xi,yi,interp=interp)
py.streamplot(xi, yi, U, V, color='k', density=1.2, minlength=.5)

py.text(-1.8, 1.7, r'$u_{i1}$', fontsize=10)
py.text(.2, .8, r'$R_1, T_1$', fontsize=10, zorder=50)
py.text(-.8, 2, r'$R_2, T_2$', fontsize=10)
py.xlabel(r'$x$', labelpad=-5)
py.ylabel(r'$y$', labelpad=-5, rotation=0)
py.xlim(-r-.01,r+.01)
py.ylim(0,r+.01)

ax = py.axes()
ax.set_aspect('equal')
ax.set_yticks([0,r])
ax.set_xticks([-r,r])
py.axis('off')

c1 = py.Circle((0,0), r, color='k', fill=False, lw=1.5)
c2 = py.Circle((d,0), 1, color='w', fill=True, clip_on=False, zorder=20)
c3 = py.Circle((d,0), 1, color='k', fill=False, lw=1.5, zorder=30)
l1 = py.Line2D([-2.1,-.4],[0,0], lw=1, c='k', clip_on=False, zorder=30)
l2 = py.Line2D([1.4,2.1],[0,0], lw=1, c='k', clip_on=False, zorder=30)
py.gca().add_artist(c1)
py.gca().add_artist(c2)
py.gca().add_artist(c3)
py.gca().add_artist(l1)
py.gca().add_artist(l2)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(py.axes())
cax = divider.append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax, ticks=[0, 1.3, 2.6])
py.text(1.2, .8, r'$\times10^{-3}$', fontsize=9)

py.tick_params(axis='both', direction='out',)
py.savefig('U.pdf', bbox_inches='tight')
