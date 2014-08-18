#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [7,3]}
py.rcParams.update(params)
import numpy as np
import sys

d = .5
r = 2

X,Y,U,V = py.loadtxt(sys.argv[1] + '.raw').T
magU = np.sqrt(U*U+V*V)
factor = 1 - int(np.log10(magU.max()))
magU *= 10**factor
print "Magnitude of U (min, max):", magU.min(), magU.max()

import matplotlib.tri as tri
triang = tri.Triangulation(X,Y)

def sqr(x):
	return x*x

xmid = X[triang.triangles].mean(axis=1)
ymid = Y[triang.triangles].mean(axis=1)
mask = np.where(sqr(xmid-d) + sqr(ymid) < sqr(1), 1, 0)
triang.set_mask(mask)

lev = np.linspace(0, (int(magU.max()*10) + 1.)/10, 11)
#lev = np.linspace(0, 2.6, 7)
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
py.streamplot(xi, yi, U, V, color='k', density=1.5, minlength=.5)

py.text(1.5, 1.7, r'$u_{i1}\times10^{' + str(factor) + r'}$')
py.text(.2, .8, r'$R_1, T_1$', zorder=50)
py.text(-1.6, 1.7, r'$R_2, T_2$')
py.xlim(-r-.01,r+.01)
py.ylim(0,r+.01)
py.axis('off')

c1 = py.Circle((0,0), r, color='k', fill=False, lw=1.5)
c2 = py.Circle((d,0), 1, color='w', fill=True, zorder=20)
c3 = py.Circle((d,0), 1, color='k', fill=False, lw=1.5, zorder=30)
l1 = py.Line2D([-2.1,-.4],[0,0], lw=1, c='k', clip_on=False, zorder=30)
l2 = py.Line2D([1.4,2.1],[0,0], lw=1, c='k', clip_on=False, zorder=30)
l3 = py.Line2D([-.4,1.4],[0,0], lw=1, c='w', clip_on=False, zorder=20)
py.gca().add_artist(c1)
py.gca().add_artist(c2)
py.gca().add_artist(c3)
py.gca().add_artist(l1)
py.gca().add_artist(l2)
py.gca().add_artist(l3)

from mpl_toolkits.axes_grid1 import make_axes_locatable
divider = make_axes_locatable(py.axes())
cax = divider.append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax)

py.tick_params(axis='both', direction='out',)
py.savefig(sys.argv[1] + '.pdf', bbox_inches='tight')
