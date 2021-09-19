#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'font.size': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np

d = 0.5
xmin, xmax = -np.pi/2, np.pi/2
def plot(name, label):
    x,y,z,s1,s2,s3 = py.loadtxt(name + 'Force.raw').T
    newdata = np.dstack((np.arctan(y/(x-d)), s1))
    xnew, ynew = np.sort(np.dstack((np.arctan2(x-d,y), s1))[0].T, axis=1);
    py.plot(xnew, np.poly1d(np.polyfit(xnew, ynew, 12))(xnew), '-', lw=1.5, label=label)
    py.plot(xnew, ynew, 'k-', lw=0.1)

plot('cylinder-inner',     r'$\mathrm{cylinders}\;\tau=4$')
plot('cylinder-outer',     r'$\mathrm{outer\;cylinder}$')
plot('sphere-inner',       r'$\mathrm{spheres}\;\tau=4$')
plot('sphere-outer',       r'$\mathrm{outer\;spheres}$')
plot('cylinder-inner-inv', r'$\mathrm{cylinders}\;\tau=-4$')
#plot('sphere-inv',   r'$\mathrm{spheres}\;\tau=-4$')

py.xlabel(r'$\varphi$', labelpad=-5)
py.ylabel(r'$-p_0 F_{x2}$', y=.75, labelpad=-5, rotation=0)
py.xlim(xmin, xmax)
ax = py.axes()
ax.set_xticks([xmin, xmax])
ax.set_xticklabels([r'$-\pi/2$', r'$\pi/2$'])
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)
ax.axhline(lw=.5, c='k', ls=':')

from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(np.pi/2))
ax.yaxis.set_minor_locator(MultipleLocator(100))

py.savefig('profiles.pdf', bbox_inches='tight', transparent=True)

