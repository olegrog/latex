#!/usr/bin/env python

import pylab as py
params = {'backend': 'pdf',
          'axes.labelsize': 10,
          'text.fontsize': 10,
          'legend.fontsize': 8,
          'xtick.labelsize': 9,
          'ytick.labelsize': 9,
          'text.usetex': True,
          'figure.figsize': [3,3]}
py.rcParams.update(params)
import numpy as np
import sys

d = float(sys.argv[2]);
xmin, xmax = -np.pi/2, np.pi/2
def plot(filename, fmt, lw, label, factor=1, correct=False):
    x,y,z,s1,s2,s3 = py.loadtxt(filename).T
    xnew, ynew = np.dstack((np.arctan2(x-d,y), s1))[0].T;
    sort = np.argsort(xnew);
    xnew, ynew = xnew[sort], ynew[sort] * factor
    if correct:
        ynew += ynew.min() * np.sin(xnew)
    py.plot(xnew, np.poly1d(np.polyfit(xnew, ynew, 12))(xnew), fmt, lw=lw, label=label)
    py.plot(xnew, ynew, 'k-', lw=0.1)

plot(sys.argv[1] + 'ShearStress1.raw', '-', 1, r'$\mathrm{term\;}p_{ij2(\gamma_1)}$')
plot(sys.argv[1] + 'ShearStress3.raw', '-', 1, r'$\mathrm{term\;}p_{ij2(\gamma_3)}$')
plot(sys.argv[1] + 'ShearStress7.raw', '-', 1, r'$\mathrm{term\;}p_{ij2(\gamma_7)}$')
plot(sys.argv[1] + 'Pressure0.raw', '-', 1, r'$\mathrm{term\;}p_2^\dag$')
plot(sys.argv[1] + 'Pressure3.raw', '-', 1, r'$\mathrm{term\;}p_{2(\gamma_3)}$')
plot(sys.argv[1] + 'Pressure7.raw', '-', 1, r'$\mathrm{term\;}p_{2(\gamma_7)}$')

py.xlabel(r'$\varphi$', labelpad=-5)
py.ylabel(r'$-F_{x2}$', y=.75, labelpad=0, rotation=0)
py.xlim(xmin, xmax)
ax = py.axes()
ax.set_xticks([xmin, xmax])
ax.set_xticklabels([r'$-\pi/2$', r'$\pi/2$'])
from matplotlib.ticker import MultipleLocator
ax.xaxis.set_minor_locator(MultipleLocator(np.pi/2))

if sys.argv[1] == 'cylinder':
    ax.yaxis.set_major_locator(MultipleLocator(200))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
if sys.argv[1] == 'sphere':
    ax.set_yticks([-300, 0, 200, 500])
    ax.yaxis.set_minor_locator(MultipleLocator(100))
if sys.argv[1] == 'cylinder-inv':
    ax.yaxis.set_major_locator(MultipleLocator(40))
    ax.yaxis.set_minor_locator(MultipleLocator(10))
if sys.argv[1] == 'sphere-inv':
    ax.yaxis.set_major_locator(MultipleLocator(40))
    ax.yaxis.set_minor_locator(MultipleLocator(100))
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc=2)


py.savefig('terms-' + sys.argv[1] + '.pdf', bbox_inches='tight')

