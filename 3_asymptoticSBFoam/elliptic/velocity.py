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

a1 = 1.5
a0 = 0.3
b0 = 0.7
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
mask = np.where(sqr(xmid/a0) + sqr(ymid/b0) < sqr(1), 1, 0)
triang.set_mask(mask)

#lev = np.linspace(magU.min(), magU.max(), 8)
lev = np.linspace(0, 12, 7)
cmap = py.cm.get_cmap('coolwarm')
CF = py.tricontourf(triang, magU, cmap=cmap) #, levels=lev)
#py.tricontour(X, Y, magU, levels=lev, colors='k', linewidths=.5)

from matplotlib.mlab import griddata
N = 50
interp = 'linear'
xi = np.linspace(0, a1, 1.5*N)
yi = np.linspace(0 , 1, N  )
U = griddata(X,Y,U,xi,yi,interp=interp)
V = griddata(X,Y,V,xi,yi,interp=interp)
py.streamplot(xi, yi, U, V, color='k', density=1.2, minlength=.5)

py.text(-.1, .9, r'$u_{i1}$', fontsize=10)
py.text(.15, .4, r'$T_1$', fontsize=10, zorder=50)
py.text(1., .8, r'$T_2$', fontsize=10)
py.xlabel(r'$x$', labelpad=-3)
py.ylabel(r'$y$', labelpad=-3, rotation=0)
py.xlim(0,a1+.05)
py.ylim(0,1+.05)

ax = py.axes()
ax.set_aspect('equal')
py.yticks([0,b0,1], (r'$0$', r'$b_0$', r'$1$'))
py.xticks([0,a0,a1], (r'$0$', r'$a_0$', r'$a_1$'))
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

from matplotlib.patches import Ellipse
c1 = Ellipse(xy=[0,0], width=2*a1, height=2, angle=0, color='k', fill=False, lw=1.5)
c2 = Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=0, color='w', fill=True, zorder=20)
c3 = Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=0, color='k', fill=False, lw=1.5, zorder=30)
py.gca().add_artist(c1)
py.gca().add_artist(c2)
py.gca().add_artist(c3)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(py.axes())
cax = divider.append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax, ticks=[0, 6, 12])
py.text(1.2, .8, r'$\times10^{-3}$', fontsize=9)

py.tick_params(axis='both', direction='out',)
py.savefig('U.pdf', bbox_inches='tight')
