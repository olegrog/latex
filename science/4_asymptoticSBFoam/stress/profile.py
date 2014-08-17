#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 10,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np
from scipy.interpolate import spline

d = 0.5
xmin, xmax = -np.pi/2, np.pi/2
def plot(filename, fmt, lw, label):
    x,y,z,s1,s2,s3 = py.loadtxt(filename).T
    newdata = np.dstack((np.arctan(y/(x-d)), s1))
    xnew, ynew = np.sort(np.dstack((np.arctan2(x-d,y), s1))[0].T, axis=1);
    py.plot(xnew, np.poly1d(np.polyfit(xnew, ynew, 12))(xnew), fmt, lw=lw, label=label)
    py.plot(xnew, ynew, 'k-', lw=0.1)

plot('cylinderShearStress.raw', 'b-', 2, r'$\mathrm{cylinders}$')
plot('sphereShearStress.raw', 'g-', 2, r'$\mathrm{spheres}$')

py.xlabel(r'$\varphi$', labelpad=-10)
py.ylabel(r'$\displaystyle\frac{p_0\tau_{ij2}d_in_j}{d}$', labelpad=-10)
py.xlim(xmin, xmax)
ax = py.axes()
ax.set_xticks([xmin,xmax])
ax.set_xticklabels([r'$-\displaystyle\frac\pi2$', r'$\displaystyle\frac\pi2$'])
ax.set_yticks([-200,0,200,400])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(np.pi/2))
ax.yaxis.set_minor_locator(MultipleLocator(100))

py.savefig('profile.pdf', bbox_inches='tight')

