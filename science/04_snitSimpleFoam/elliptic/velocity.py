#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'font.size': 10,
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [7.5,4]}
py.rcParams.update(params)
import numpy as np
from matplotlib.tri import Triangulation
from matplotlib.mlab import griddata
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Ellipse
from matplotlib.patches import Arc

a1 = 1.5
a0 = 0.3
b0 = 0.7
beta = np.pi/4

def draw(filename, ax, beta, lw):
    X,Y,U,V = py.loadtxt(filename).T
    magU = np.sqrt(U*U+V*V)

    triang = Triangulation(X,Y)
    xmid = X[triang.triangles].mean(axis=1)
    ymid = Y[triang.triangles].mean(axis=1)

    lev = np.linspace(0, (int(magU.max()*10) + 1.)/10, 11)
    cmap = py.cm.get_cmap('coolwarm')
    CF = ax.tricontourf(triang, magU, cmap=cmap, levels=lev)

    N = 50
    interp = 'linear'
    xi = np.linspace(X.min(), X.max(), (X.max()-X.min())*N)
    yi = np.linspace(Y.min(), Y.max(), (Y.max()-Y.min())*N)
    U = griddata(X,Y,U,xi,yi,interp=interp)
    V = griddata(X,Y,V,xi,yi,interp=interp)
    ax.streamplot(xi, yi, U, V, color='k', density=2*np.sqrt(Y.max()-Y.min()), minlength=.2, linewidth=lw, arrowstyle='->', arrowsize=lw)

    ax.set_aspect('equal')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()

    angle = -(np.pi/2-beta)/np.pi*180
    artists = [
        Ellipse(xy=[0,0], width=2*a1, height=2, angle=0, color='k', fill=False, lw=1.5*lw),
        Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=angle, color='w', fill=True, zorder=20),
        Ellipse(xy=[0,0], width=2*a0, height=2*b0, angle=angle, color='k', fill=False, lw=1.5*lw, zorder=30),
    ]

    for a in artists:
        ax.add_artist(a)
    
    return CF

gs = GridSpec(1, 2, width_ratios=[2.2, 1], wspace=0.25)
ax1 = py.subplot(gs[0])
ax2 = py.subplot(gs[1])
CF = draw('elliptic.raw', ax1, np.pi/2, 1)
draw('rotated.raw', ax2, beta, .7)
cax = make_axes_locatable(ax1).append_axes("right", "5%", pad="5%")
py.colorbar(CF, cax=cax)

### Subplot 1 preferences
py.setp(ax1,
    xticks=[0, a0, a1], xticklabels=[r'$0$', r'$a_0$', r'$a_1$'],
    yticks=[0, b0, 1], yticklabels=[r'$0$', r'$b_0$', r'$1$'])

ax1.text(1.5, .9, r'$u_{i1}$')
ax1.text(.15, .4, r'$T_1$', zorder=50)
ax1.text(1., .8, r'$T_2$')
ax1.set_xlim(0, a1+.05)
ax1.set_ylim(0, 1+.05)
ax1.text(a1/2, -.09, r'$x$')
ax1.text(-.09, .85, r'$y$')
ax1.set_xlabel(r'$\mathrm{(a)}$', labelpad=10)

### Subplot 2 preferences
py.setp(ax2,
    xticks=[0, a1], xticklabels=[r'$0$', r'$a_1$'],
    yticks=[-1, 0, 1], yticklabels=[r'$-1$', r'$0$', r'$1$'])

ax2.set_xlim(0, a1+.01)
ax2.set_ylim(-1, 1.035)
ax2.add_artist(py.Line2D([-.03, b0/np.sqrt(2)+.03],[-.03, b0/np.sqrt(2)+.03], lw=.5, c='k', clip_on=False, zorder=30))
ax2.add_artist(py.Line2D([-.03, .42],[0, 0], lw=.5, c='k', clip_on=False, zorder=30))
ax2.add_artist(Arc(xy=[0,0], width=.3, height=.3, theta1=0, theta2=45, color='k', zorder=30))
ax2.text(.2, .07, r'$\beta$', zorder=50)
ax2.text(a1/2, -1.17, r'$x$')
ax2.text(-0.17, .75, r'$y$')
ax2.set_xlabel(r'$\mathrm{(b)}$', labelpad=10)

py.tick_params(axis='both', direction='out',)
py.savefig('U.pdf', bbox_inches='tight', transparent=True)
